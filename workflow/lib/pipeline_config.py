#!/usr/bin/env python3
"""Configuration and sample-sheet validation shared by both workflow engines."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Any

import yaml

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,99}$")
MODES = {"bulk_rna", "dna_reseq", "denovo_genome", "denovo_transcriptome"}
LAYOUTS = {"paired", "single"}
STRANDEDNESS = {"unstranded", "forward", "reverse"}
PLACEHOLDER_TOKENS = ("/GROUP/", "YOUR_", "REPLACE_", "<", ">")
STAR_INDEX_FILES = ("Genome", "SA", "SAindex", "chrLength.txt", "chrName.txt")
BWA_MEM2_INDEX_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")


class ConfigError(ValueError):
    """Raised when a project configuration is incomplete or unsafe."""


def _require(mapping: dict[str, Any], key: str, context: str) -> Any:
    value = mapping.get(key)
    if value is None or value == "":
        raise ConfigError(f"Missing required setting: {context}.{key}")
    return value


def _resolve_path(value: str | Path) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(str(value)))
    return Path(expanded).resolve(strict=False)


def _check_placeholder(path: Path, label: str) -> None:
    value = str(path)
    if any(token in value for token in PLACEHOLDER_TOKENS):
        raise ConfigError(f"{label} still contains an example placeholder: {value}")


def _check_file(path: Path, label: str, check_files: bool) -> None:
    _check_placeholder(path, label)
    if check_files and not path.is_file():
        raise ConfigError(f"{label} is not a readable file: {path}")
    if check_files and not os.access(path, os.R_OK):
        raise ConfigError(f"{label} is not readable: {path}")


def _read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ConfigError(f"TSV has no header: {path}")
        rows = []
        for row in reader:
            if not any((value or "").strip() for value in row.values()):
                continue
            rows.append({key: (value or "").strip() for key, value in row.items()})
    return rows, list(reader.fieldnames)


def _load_samples(
    sample_sheet: Path,
    layout: str,
    mode: str,
    design: list[str],
    check_files: bool,
) -> tuple[dict[str, dict[str, str]], list[str]]:
    _check_file(sample_sheet, "project.sample_sheet", check_files)
    rows, columns = _read_tsv(sample_sheet)
    required = {"sample_id", "read1"}
    if layout == "paired":
        required.add("read2")
    if mode == "bulk_rna":
        required.update(design)
    missing = sorted(required.difference(columns))
    if missing:
        raise ConfigError(f"Sample sheet is missing columns: {', '.join(missing)}")
    if not rows:
        raise ConfigError("Sample sheet contains no samples")

    samples: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        sample_id = row.get("sample_id", "")
        if not SAFE_ID.fullmatch(sample_id):
            raise ConfigError(
                f"Invalid sample_id on line {line_number}: {sample_id!r}. "
                "Use letters, digits, dot, underscore, or hyphen."
            )
        if sample_id in samples:
            raise ConfigError(f"Duplicate sample_id: {sample_id}")

        read1 = _resolve_path(row.get("read1", ""))
        _check_file(read1, f"{sample_id}.read1", check_files)
        row["read1"] = str(read1)

        read2_value = row.get("read2", "")
        if layout == "paired":
            if not read2_value:
                raise ConfigError(f"Paired sample {sample_id} has no read2 path")
            read2 = _resolve_path(read2_value)
            _check_file(read2, f"{sample_id}.read2", check_files)
            if read1 == read2:
                raise ConfigError(f"Sample {sample_id} uses the same file for read1 and read2")
            row["read2"] = str(read2)
        else:
            if read2_value:
                raise ConfigError(
                    f"Single-end run has a read2 value for {sample_id}; use one layout per run"
                )
            row["read2"] = ""

        for term in design:
            if not row.get(term, ""):
                raise ConfigError(f"Sample {sample_id} has an empty design value for {term}")
        samples[sample_id] = row

    return samples, columns


def _load_contrasts(
    contrast_path: Path | None,
    columns: list[str],
    samples: dict[str, dict[str, str]],
    check_files: bool,
) -> list[dict[str, str]]:
    if contrast_path is None:
        return []
    _check_file(contrast_path, "rna.contrasts", check_files)
    rows, contrast_columns = _read_tsv(contrast_path)
    required = {"contrast_id", "factor", "numerator", "denominator"}
    missing = sorted(required.difference(contrast_columns))
    if missing:
        raise ConfigError(f"Contrast table is missing columns: {', '.join(missing)}")

    contrasts: list[dict[str, str]] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        contrast_id = row["contrast_id"]
        if not SAFE_ID.fullmatch(contrast_id):
            raise ConfigError(f"Invalid contrast_id on line {line_number}: {contrast_id!r}")
        if contrast_id in seen:
            raise ConfigError(f"Duplicate contrast_id: {contrast_id}")
        seen.add(contrast_id)

        factor = row["factor"]
        if factor not in columns:
            raise ConfigError(f"Contrast factor {factor!r} is not a sample-sheet column")
        values = {sample[factor] for sample in samples.values()}
        if row["numerator"] not in values:
            raise ConfigError(
                f"Contrast {contrast_id} numerator {row['numerator']!r} is absent from {factor}"
            )
        if row["denominator"] not in values:
            raise ConfigError(
                f"Contrast {contrast_id} denominator {row['denominator']!r} is absent from {factor}"
            )
        if row["numerator"] == row["denominator"]:
            raise ConfigError(f"Contrast {contrast_id} compares a group with itself")
        if not row.get("label"):
            row["label"] = contrast_id
        contrasts.append(row)
    return contrasts


def validate_and_load(config: dict[str, Any], check_files: bool = True) -> dict[str, Any]:
    """Validate configuration and return normalized workflow state."""

    project = config.get("project", {})
    analysis = config.get("analysis", {})
    reference = config.get("reference", {})
    qc = config.get("qc", {})
    rna = config.get("rna", {})
    assembly = config.get("assembly", {})

    project_id = str(_require(project, "id", "project"))
    if not SAFE_ID.fullmatch(project_id):
        raise ConfigError(
            "project.id must contain only letters, digits, dot, underscore, or hyphen"
        )

    output_dir = _resolve_path(_require(project, "output_dir", "project"))
    _check_placeholder(output_dir, "project.output_dir")
    if output_dir == Path("/"):
        raise ConfigError("project.output_dir cannot be the filesystem root")
    if check_files:
        if output_dir.exists() and not output_dir.is_dir():
            raise ConfigError(f"project.output_dir exists but is not a directory: {output_dir}")
        writable_parent = output_dir
        while not writable_parent.exists() and writable_parent != writable_parent.parent:
            writable_parent = writable_parent.parent
        if not os.access(writable_parent, os.W_OK | os.X_OK):
            raise ConfigError(
                f"project.output_dir is not creatable under {writable_parent}: {output_dir}"
            )

    sample_sheet = _resolve_path(_require(project, "sample_sheet", "project"))
    mode = str(_require(analysis, "mode", "analysis"))
    layout = str(_require(analysis, "layout", "analysis"))
    if mode not in MODES:
        raise ConfigError(f"analysis.mode must be one of: {', '.join(sorted(MODES))}")
    if layout not in LAYOUTS:
        raise ConfigError(f"analysis.layout must be one of: {', '.join(sorted(LAYOUTS))}")

    design = list(rna.get("design", ["condition"])) if mode == "bulk_rna" else []
    if mode == "bulk_rna":
        valid_terms = all(
            isinstance(term, str) and SAFE_ID.fullmatch(term) for term in design
        )
        if not design or not valid_terms:
            raise ConfigError(
                "rna.design must be a non-empty list of safe sample-sheet column names"
            )
        strandedness = str(rna.get("strandedness", "unstranded"))
        if strandedness not in STRANDEDNESS:
            raise ConfigError(
                f"rna.strandedness must be one of: {', '.join(sorted(STRANDEDNESS))}"
            )

    samples, sample_columns = _load_samples(
        sample_sheet, layout, mode, design, check_files
    )
    if mode == "bulk_rna":
        for term in design:
            levels = {sample[term] for sample in samples.values()}
            if len(levels) < 2:
                raise ConfigError(
                    f"RNA design term {term!r} has fewer than two observed levels"
                )

    fasta: Path | None = None
    gtf: Path | None = None
    if mode in {"bulk_rna", "dna_reseq"}:
        fasta = _resolve_path(_require(reference, "fasta", "reference"))
        _check_file(fasta, "reference.fasta", check_files)
    if mode == "bulk_rna":
        gtf = _resolve_path(_require(reference, "annotation_gtf", "reference"))
        _check_file(gtf, "reference.annotation_gtf", check_files)

    star_index: Path | None = None
    if mode == "bulk_rna" and reference.get("star_index"):
        star_index = _resolve_path(reference["star_index"])
        _check_placeholder(star_index, "reference.star_index")
        missing = [name for name in STAR_INDEX_FILES if not (star_index / name).is_file()]
        if check_files and missing:
            raise ConfigError(
                f"STAR index is incomplete at {star_index}; missing: {', '.join(missing)}"
            )

    bwa_prefix: Path | None = None
    if mode == "dna_reseq" and reference.get("bwa_index_prefix"):
        bwa_prefix = _resolve_path(reference["bwa_index_prefix"])
        _check_placeholder(bwa_prefix, "reference.bwa_index_prefix")
        missing = [
            f"{bwa_prefix}{suffix}"
            for suffix in BWA_MEM2_INDEX_SUFFIXES
            if not Path(f"{bwa_prefix}{suffix}").is_file()
        ]
        if check_files and missing:
            raise ConfigError(
                f"BWA-MEM2 index is incomplete; missing: {', '.join(missing)}"
            )

    contrast_path: Path | None = None
    contrasts: list[dict[str, str]] = []
    if mode == "bulk_rna" and rna.get("contrasts"):
        contrast_path = _resolve_path(rna["contrasts"])
        contrasts = _load_contrasts(
            contrast_path, sample_columns, samples, check_files
        )
        for contrast in contrasts:
            if contrast["factor"] not in design:
                raise ConfigError(
                    f"Contrast {contrast['contrast_id']} factor {contrast['factor']!r} "
                    "is not listed in rna.design"
                )
            numerator_n = sum(
                sample[contrast["factor"]] == contrast["numerator"]
                for sample in samples.values()
            )
            denominator_n = sum(
                sample[contrast["factor"]] == contrast["denominator"]
                for sample in samples.values()
            )
            if numerator_n < 2 or denominator_n < 2:
                raise ConfigError(
                    f"Contrast {contrast['contrast_id']} requires at least two samples "
                    "in both numerator and denominator groups"
                )

    screen = qc.get("fastq_screen", {})
    screen_config: Path | None = None
    if bool(screen.get("enabled", False)):
        screen_config = _resolve_path(_require(screen, "config", "qc.fastq_screen"))
        _check_file(screen_config, "qc.fastq_screen.config", check_files)

    busco = assembly.get("busco", {})
    busco_lineage: Path | None = None
    if mode.startswith("denovo_") and bool(busco.get("enabled", False)):
        busco_lineage = _resolve_path(_require(busco, "lineage_path", "assembly.busco"))
        _check_placeholder(busco_lineage, "assembly.busco.lineage_path")
        if check_files and not busco_lineage.is_dir():
            raise ConfigError(f"BUSCO lineage directory does not exist: {busco_lineage}")

    return {
        "project_id": project_id,
        "display_name": str(project.get("display_name", project_id)),
        "output_dir": str(output_dir),
        "sample_sheet": str(sample_sheet),
        "samples": samples,
        "sample_columns": sample_columns,
        "sample_ids": list(samples),
        "mode": mode,
        "layout": layout,
        "paired": layout == "paired",
        "fasta": str(fasta) if fasta else None,
        "gtf": str(gtf) if gtf else None,
        "star_index": str(star_index) if star_index else None,
        "bwa_index_prefix": str(bwa_prefix) if bwa_prefix else None,
        "design": design,
        "contrast_path": str(contrast_path) if contrast_path else None,
        "contrasts": contrasts,
        "contrast_ids": [row["contrast_id"] for row in contrasts],
        "screen_enabled": bool(screen.get("enabled", False)),
        "screen_config": str(screen_config) if screen_config else None,
        "busco_enabled": bool(busco.get("enabled", False)),
        "busco_lineage": str(busco_lineage) if busco_lineage else None,
    }


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--no-check-files", action="store_true")
    args = parser.parse_args()
    try:
        with args.config.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle) or {}
        state = validate_and_load(loaded, check_files=not args.no_check_files)
    except (OSError, yaml.YAMLError, ConfigError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print("Configuration valid")
    print(f"  project: {state['project_id']}")
    print(f"  mode: {state['mode']}")
    print(f"  layout: {state['layout']}")
    print(f"  samples: {len(state['sample_ids'])}")
    print(f"  output: {state['output_dir']}")
    if state["contrast_ids"]:
        print(f"  contrasts: {', '.join(state['contrast_ids'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
