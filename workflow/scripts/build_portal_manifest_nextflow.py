#!/usr/bin/env python3
"""Build the Shiny manifest from staged Nextflow outputs without exposing inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def percent(value: Any) -> float | None:
    return None if value is None else round(float(value) * 100, 4)


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
        "q30_percent": percent(after.get("q30_rate")),
        "gc_percent": percent(after.get("gc_content")),
        "duplication_percent": percent(duplication.get("rate")),
        "adapter_trimmed_reads": adapter.get("adapter_trimmed_reads"),
    }


def parse_star(path: Path, sample_id: str) -> dict[str, Any]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "|" in line:
            key, value = (part.strip() for part in line.split("|", 1))
            values[key] = value

    def number(key: str) -> float | None:
        raw = values.get(key)
        if raw is None:
            return None
        try:
            return float(raw.replace("%", "").replace(",", ""))
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
    row: dict[str, Any] = {
        "sample_id": sample_id,
        "input_reads": None,
        "mapped_reads": None,
        "mapping_percent": None,
        "multimapping_percent": None,
        "unmapped_too_short_percent": None,
        "mean_coverage": None,
    }
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if " in total " in line:
            row["input_reads"] = float(line.split("+", 1)[0].strip())
        elif " mapped (" in line and "primary mapped" not in line:
            row["mapped_reads"] = float(line.split("+", 1)[0].strip())
            try:
                row["mapping_percent"] = float(line.split("(", 1)[1].split("%", 1)[0])
            except (IndexError, ValueError):
                pass
    return row


def parse_mosdepth(path: Path) -> float | None:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row.get("chrom") in {"total", "total_region"}:
            try:
                return float(row["mean"])
            except (KeyError, TypeError, ValueError):
                return None
    return None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build(args: argparse.Namespace) -> None:
    config = json.loads(args.config.read_text(encoding="utf-8"))
    samples = pd.read_csv(args.sample_sheet, sep="\t", dtype=str).fillna("")
    sample_ids = samples["sample_id"].tolist()
    output_root = args.output_root
    portal = output_root / "portal"
    tables = portal / "tables"
    provenance = getattr(args, "provenance_dir", None) or output_root / "provenance"
    tables.mkdir(parents=True, exist_ok=True)
    provenance.mkdir(parents=True, exist_ok=True)

    design_terms = config.get("rna", {}).get("design", []) if args.mode == "bulk_rna" else []
    visible = list(dict.fromkeys(["sample_id", *[x for x in design_terms if x in samples.columns]]))
    samples_path = tables / "samples.tsv"
    samples[visible].to_csv(samples_path, sep="\t", index=False)

    qc_rows = []
    for sample_id in sample_ids:
        path = args.fastp_dir / f"{sample_id}.fastp.json"
        if path.is_file():
            qc_rows.append(parse_fastp(path, sample_id))
    qc_path = tables / "sample_qc.tsv"
    pd.DataFrame(qc_rows or [{"sample_id": sample_ids[0]}]).iloc[: len(qc_rows)].to_csv(
        qc_path, sep="\t", index=False
    )

    mapping_rows = []
    if args.mode == "bulk_rna":
        for sample_id in sample_ids:
            path = args.analysis_dir / f"{sample_id}.Log.final.out"
            if path.is_file():
                mapping_rows.append(parse_star(path, sample_id))
    elif args.mode == "dna_reseq":
        for sample_id in sample_ids:
            flagstat = args.analysis_dir / f"{sample_id}.flagstat.txt"
            if flagstat.is_file():
                row = parse_flagstat(flagstat, sample_id)
                mosdepth = args.analysis_dir / f"{sample_id}.mosdepth.summary.txt"
                if mosdepth.is_file():
                    row["mean_coverage"] = parse_mosdepth(mosdepth)
                mapping_rows.append(row)
    mapping_path = tables / "mapping_qc.tsv"
    pd.DataFrame(mapping_rows or [{"sample_id": sample_ids[0]}]).iloc[: len(mapping_rows)].to_csv(
        mapping_path, sep="\t", index=False
    )

    downloads: list[dict[str, str]] = []
    checksum_sources: dict[str, Path] = {}

    def add_download(label: str, relative: str, source: Path, kind: str) -> None:
        if source.is_file():
            downloads.append({"label": label, "path": relative, "kind": kind})
            checksum_sources[relative] = source

    add_download(
        "MultiQC report",
        "results/multiqc/multiqc_report.html",
        args.multiqc,
        "html",
    )
    figures: dict[str, str] = {}
    differential_expression: list[dict[str, str]] = []
    assemblies: list[dict[str, str]] = []

    if args.mode == "bulk_rna":
        counts = args.analysis_dir / "gene_counts.tsv"
        add_download("Raw gene counts", "results/rna/counts/gene_counts.tsv", counts, "table")
        for sample_id in sample_ids:
            track = args.analysis_dir / f"{sample_id}.CPM.bw"
            add_download(
                f"{sample_id}: normalized coverage track",
                f"results/rna/tracks/{sample_id}.CPM.bw",
                track,
                "bigwig",
            )
        de_root = args.analysis_dir / "deseq2"
        add_download(
            "Normalized gene counts",
            "results/rna/deseq2/normalized_counts.tsv",
            de_root / "normalized_counts.tsv",
            "table",
        )
        add_download(
            "Contrast summary",
            "results/rna/deseq2/contrast_summary.tsv",
            de_root / "contrast_summary.tsv",
            "table",
        )
        if (de_root / "pca.png").is_file():
            figures["pca"] = "results/rna/deseq2/pca.png"
        if (de_root / "sample_distance.png").is_file():
            figures["sample_distance"] = "results/rna/deseq2/sample_distance.png"
        contrast_file = args.analysis_dir / "contrasts.tsv"
        labels: dict[str, str] = {}
        if contrast_file.is_file():
            frame = pd.read_csv(contrast_file, sep="\t", dtype=str).fillna("")
            labels = {
                row.contrast_id: row.label or row.contrast_id
                for row in frame.itertuples(index=False)
            }
        for contrast_id, label in labels.items():
            result = de_root / contrast_id / "results.tsv"
            volcano = de_root / contrast_id / "volcano.png"
            ma_plot = de_root / contrast_id / "ma.png"
            if result.is_file():
                relative_root = f"results/rna/deseq2/{contrast_id}"
                add_download(
                    f"{label}: differential expression",
                    f"{relative_root}/results.tsv",
                    result,
                    "table",
                )
                differential_expression.append(
                    {
                        "id": contrast_id,
                        "label": label,
                        "table": f"{relative_root}/results.tsv",
                        "volcano": f"{relative_root}/volcano.png" if volcano.is_file() else "",
                        "ma_plot": f"{relative_root}/ma.png" if ma_plot.is_file() else "",
                    }
                )
    elif args.mode == "dna_reseq":
        for sample_id in sample_ids:
            vcf = args.analysis_dir / f"{sample_id}.filtered.vcf.gz"
            stats = args.analysis_dir / f"{sample_id}.bcftools.stats.txt"
            add_download(
                f"{sample_id}: filtered VCF",
                f"results/dna/variants/{sample_id}.filtered.vcf.gz",
                vcf,
                "vcf",
            )
            add_download(
                f"{sample_id}: variant statistics",
                f"results/dna/variants/{sample_id}.bcftools.stats.txt",
                stats,
                "text",
            )
    else:
        for sample_id in sample_ids:
            contigs = args.analysis_dir / f"{sample_id}.contigs.fasta"
            quast_html = args.analysis_dir / f"{sample_id}.quast.report.html"
            quast_table = args.analysis_dir / f"{sample_id}.quast.report.tsv"
            base = f"results/assembly/{sample_id}"
            add_download(
                f"{sample_id}: assembly FASTA",
                f"{base}/{sample_id}.contigs.fasta",
                contigs,
                "fasta",
            )
            add_download(
                f"{sample_id}: QUAST report",
                f"{base}/quast/{sample_id}.quast.report.html",
                quast_html,
                "html",
            )
            assemblies.append(
                {
                    "sample_id": sample_id,
                    "contigs": f"{base}/{sample_id}.contigs.fasta" if contigs.is_file() else "",
                    "quast_html": (
                        f"{base}/quast/{sample_id}.quast.report.html"
                        if quast_html.is_file()
                        else ""
                    ),
                    "quast_table": (
                        f"{base}/quast/{sample_id}.quast.report.tsv"
                        if quast_table.is_file()
                        else ""
                    ),
                }
            )

    manifest = {
        "schema_version": "1.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "pipeline": {
            "name": "Hellbender NGS Pipeline",
            "version": str(config.get("pipeline", {}).get("version", "1.1.0")),
            "execution_engine": "Nextflow",
        },
        "project": {
            "id": str(config["project"]["id"]),
            "display_name": str(config["project"].get("display_name", config["project"]["id"])),
            "analysis_mode": args.mode,
            "description": str(config.get("portal", {}).get("description", "")),
            "status": "complete",
        },
        "tables": {
            "samples": "portal/tables/samples.tsv",
            "sample_qc": "portal/tables/sample_qc.tsv",
            "mapping_qc": "portal/tables/mapping_qc.tsv",
        },
        "figures": figures,
        "differential_expression": differential_expression,
        "assemblies": assemblies,
        "downloads": downloads,
    }
    (portal / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (provenance / "config.used.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )

    checksum_sources.update(
        {
            "portal/tables/samples.tsv": samples_path,
            "portal/tables/sample_qc.tsv": qc_path,
            "portal/tables/mapping_qc.tsv": mapping_path,
        }
    )
    lines = []
    for relative, source in sorted(checksum_sources.items()):
        if source.is_file() and source.stat().st_size <= 250 * 1024 * 1024:
            lines.append(f"{sha256(source)}  {relative}")
    (provenance / "result_checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--fastp-dir", required=True, type=Path)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--multiqc", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--provenance-dir", type=Path)
    args = parser.parse_args()
    build(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
