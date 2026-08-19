"""Build a privacy-conscious, client-facing manifest from completed outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

project_dir = Path(str(snakemake.params.project_dir)).resolve()  # noqa: F821
project_id = str(snakemake.params.project_id)  # noqa: F821
display_name = str(snakemake.params.display_name)  # noqa: F821
mode = str(snakemake.params.mode)  # noqa: F821
pipeline_version = str(snakemake.params.pipeline_version)  # noqa: F821
sample_sheet = Path(str(snakemake.params.sample_sheet)).resolve()  # noqa: F821
contrast_ids = list(snakemake.params.contrast_ids)  # noqa: F821
contrast_rows = list(snakemake.params.contrasts)  # noqa: F821
config = dict(snakemake.params.config)  # noqa: F821

manifest_path = Path(str(snakemake.output.manifest))  # noqa: F821
samples_path = Path(str(snakemake.output.samples))  # noqa: F821
qc_path = Path(str(snakemake.output.qc))  # noqa: F821
mapping_path = Path(str(snakemake.output.mapping))  # noqa: F821
config_used_path = Path(str(snakemake.output.config_used))  # noqa: F821
checksums_path = Path(str(snakemake.output.checksums))  # noqa: F821
log_path = Path(str(snakemake.log[0]))  # noqa: F821

for path in (
    manifest_path,
    samples_path,
    qc_path,
    mapping_path,
    config_used_path,
    checksums_path,
    log_path,
):
    path.parent.mkdir(parents=True, exist_ok=True)


def relpath(path: Path) -> str:
    """Return a POSIX relative path only when it stays inside the project."""

    resolved = path.resolve(strict=False)
    try:
        return resolved.relative_to(project_dir).as_posix()
    except ValueError as exc:
        raise ValueError(f"Portal file is outside the project directory: {resolved}") from exc


def parse_fastp(path: Path, sample_id: str) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    duplication = data.get("duplication", {})
    adapter = data.get("adapter_cutting", {})
    return {
        "sample_id": sample_id,
        "reads_before": before.get("total_reads"),
        "reads_after": after.get("total_reads"),
        "bases_after": after.get("total_bases"),
        "q30_percent": _percent(after.get("q30_rate")),
        "gc_percent": _percent(after.get("gc_content")),
        "duplication_percent": _percent(duplication.get("rate")),
        "adapter_trimmed_reads": adapter.get("adapter_trimmed_reads"),
    }


def _percent(value: Any) -> float | None:
    if value is None:
        return None
    return round(float(value) * 100, 4)


def parse_star(path: Path, sample_id: str) -> dict[str, Any]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "|" not in line:
            continue
        key, value = (part.strip() for part in line.split("|", 1))
        values[key] = value

    def number(key: str) -> float | None:
        raw = values.get(key)
        if raw is None:
            return None
        cleaned = raw.replace("%", "").replace(",", "")
        try:
            return float(cleaned)
        except ValueError:
            return None

    return {
        "sample_id": sample_id,
        "input_reads": number("Number of input reads"),
        "mapped_reads": number("Uniquely mapped reads number"),
        "mapping_percent": number("Uniquely mapped reads %"),
        "multimapping_percent": number("% of reads mapped to multiple loci"),
        "unmapped_too_short_percent": number("% of reads unmapped: too short"),
        "mean_coverage": None,
    }


def parse_flagstat(path: Path, sample_id: str) -> dict[str, Any]:
    metrics: dict[str, float | None] = {
        "input_reads": None,
        "mapped_reads": None,
        "mapping_percent": None,
        "multimapping_percent": None,
        "unmapped_too_short_percent": None,
        "mean_coverage": None,
    }
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if " in total " in line:
            metrics["input_reads"] = float(line.split("+", 1)[0].strip())
        elif " mapped (" in line and "primary mapped" not in line:
            metrics["mapped_reads"] = float(line.split("+", 1)[0].strip())
            if "(" in line and "%" in line:
                try:
                    metrics["mapping_percent"] = float(line.split("(", 1)[1].split("%", 1)[0])
                except ValueError:
                    pass
    return {"sample_id": sample_id, **metrics}


def parse_mosdepth_mean(path: Path) -> float | None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    for row in rows:
        if row.get("chrom") == "total_region" or row.get("chrom") == "total":
            try:
                return float(row["mean"])
            except (KeyError, TypeError, ValueError):
                return None
    return None


samples = pd.read_csv(sample_sheet, sep="\t", dtype=str).fillna("")
design_terms = list(config.get("rna", {}).get("design", [])) if mode == "bulk_rna" else []
visible_columns = ["sample_id", *[term for term in design_terms if term in samples.columns]]
visible_columns = list(dict.fromkeys(visible_columns))
samples[visible_columns].to_csv(samples_path, sep="\t", index=False)

qc_rows = []
for sample_id in samples["sample_id"]:
    fastp_path = project_dir / "results" / "qc" / "fastp" / f"{sample_id}.fastp.json"
    if fastp_path.is_file():
        qc_rows.append(parse_fastp(fastp_path, sample_id))
qc_frame = pd.DataFrame(qc_rows)
if qc_frame.empty:
    qc_frame = pd.DataFrame(columns=["sample_id"])
qc_frame.to_csv(qc_path, sep="\t", index=False)

mapping_rows = []
for sample_id in samples["sample_id"]:
    if mode == "bulk_rna":
        star_log = project_dir / "results" / "rna" / "alignment" / f"{sample_id}.Log.final.out"
        if star_log.is_file():
            mapping_rows.append(parse_star(star_log, sample_id))
    elif mode == "dna_reseq":
        flagstat = project_dir / "results" / "dna" / "qc" / f"{sample_id}.flagstat.txt"
        if flagstat.is_file():
            row = parse_flagstat(flagstat, sample_id)
            mosdepth = project_dir / "results" / "dna" / "qc" / f"{sample_id}.mosdepth.summary.txt"
            if mosdepth.is_file():
                row["mean_coverage"] = parse_mosdepth_mean(mosdepth)
            mapping_rows.append(row)
mapping_frame = pd.DataFrame(mapping_rows)
if mapping_frame.empty:
    mapping_frame = pd.DataFrame(columns=["sample_id"])
mapping_frame.to_csv(mapping_path, sep="\t", index=False)

downloads: list[dict[str, str]] = []
figures: dict[str, Any] = {}
differential_expression: list[dict[str, str]] = []
assemblies: list[dict[str, str]] = []


def add_download(label: str, path: Path, kind: str) -> None:
    if path.is_file():
        downloads.append({"label": label, "path": relpath(path), "kind": kind})


multiqc = project_dir / "results" / "multiqc" / "multiqc_report.html"
add_download("MultiQC report", multiqc, "html")

if mode == "bulk_rna":
    counts = project_dir / "results" / "rna" / "counts" / "gene_counts.tsv"
    add_download("Raw gene counts", counts, "table")
    de_root = project_dir / "results" / "rna" / "deseq2"
    normalized = de_root / "normalized_counts.tsv"
    contrast_summary = de_root / "contrast_summary.tsv"
    add_download("Normalized gene counts", normalized, "table")
    add_download("Contrast summary", contrast_summary, "table")
    if (de_root / "pca.png").is_file():
        figures["pca"] = relpath(de_root / "pca.png")
    if (de_root / "sample_distance.png").is_file():
        figures["sample_distance"] = relpath(de_root / "sample_distance.png")
    contrast_lookup = {
        row.get("contrast_id"): row for row in contrast_rows if isinstance(row, dict)
    }
    for contrast_id in contrast_ids:
        contrast_dir = de_root / contrast_id
        result_table = contrast_dir / "results.tsv"
        volcano = contrast_dir / "volcano.png"
        ma_plot = contrast_dir / "ma.png"
        if result_table.is_file():
            add_download(f"{contrast_id}: differential expression", result_table, "table")
            differential_expression.append(
                {
                    "id": contrast_id,
                    "label": contrast_lookup.get(contrast_id, {}).get("label", contrast_id),
                    "table": relpath(result_table),
                    "volcano": relpath(volcano) if volcano.is_file() else "",
                    "ma_plot": relpath(ma_plot) if ma_plot.is_file() else "",
                }
            )

elif mode == "dna_reseq":
    for sample_id in samples["sample_id"]:
        vcf = project_dir / "results" / "dna" / "variants" / f"{sample_id}.filtered.vcf.gz"
        stats = project_dir / "results" / "dna" / "variants" / f"{sample_id}.bcftools.stats.txt"
        add_download(f"{sample_id}: filtered VCF", vcf, "vcf")
        add_download(f"{sample_id}: variant statistics", stats, "text")

else:
    for sample_id in samples["sample_id"]:
        contigs = project_dir / "results" / "assembly" / sample_id / "contigs.fasta"
        quast_html = project_dir / "results" / "assembly" / sample_id / "quast" / "report.html"
        quast_tsv = project_dir / "results" / "assembly" / sample_id / "quast" / "report.tsv"
        add_download(f"{sample_id}: assembly FASTA", contigs, "fasta")
        add_download(f"{sample_id}: QUAST report", quast_html, "html")
        assemblies.append(
            {
                "sample_id": sample_id,
                "contigs": relpath(contigs) if contigs.is_file() else "",
                "quast_html": relpath(quast_html) if quast_html.is_file() else "",
                "quast_table": relpath(quast_tsv) if quast_tsv.is_file() else "",
            }
        )

manifest = {
    "schema_version": "1.0",
    "generated_at": datetime.now(UTC).isoformat(),
    "pipeline": {"name": "Hellbender NGS Pipeline", "version": pipeline_version},
    "project": {
        "id": project_id,
        "display_name": display_name,
        "analysis_mode": mode,
        "description": str(config.get("portal", {}).get("description", "")),
        "status": "complete",
    },
    "tables": {
        "samples": relpath(samples_path),
        "sample_qc": relpath(qc_path),
        "mapping_qc": relpath(mapping_path),
    },
    "figures": figures,
    "differential_expression": differential_expression,
    "assemblies": assemblies,
    "downloads": downloads,
}

manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
config_used_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

hash_candidates = [samples_path, qc_path, mapping_path, multiqc]
hash_candidates.extend(project_dir / item["path"] for item in downloads)
seen_paths: set[Path] = set()
checksum_lines = []
for path in hash_candidates:
    path = path.resolve(strict=False)
    if path in seen_paths or not path.is_file():
        continue
    seen_paths.add(path)
    if path.stat().st_size > 250 * 1024 * 1024:
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    checksum_lines.append(f"{digest.hexdigest()}  {relpath(path)}")
checksums_path.write_text("\n".join(sorted(checksum_lines)) + "\n", encoding="utf-8")

log_path.write_text(
    "\n".join(
        [
            f"project={project_id}",
            f"mode={mode}",
            f"samples={len(samples)}",
            f"downloads={len(downloads)}",
            f"manifest={manifest_path}",
        ]
    )
    + "\n",
    encoding="utf-8",
)
