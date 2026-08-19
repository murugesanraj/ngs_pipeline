from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

from aggregate_featurecounts_nextflow import clean_featurecounts  # noqa: E402
from build_portal_manifest_nextflow import build  # noqa: E402


def test_nextflow_featurecounts_columns_are_normalized(tmp_path: Path) -> None:
    raw = tmp_path / "featurecounts.raw.tsv"
    raw.write_text(
        "# Program:featureCounts\n"
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\tbams/control_1.bam\tbams/treated_1.bam\n"
        "gene_a\tchr1\t1\t20\t+\t20\t7\t11\n",
        encoding="utf-8",
    )
    output = tmp_path / "gene_counts.tsv"

    clean_featurecounts(raw, output)

    counts = pd.read_csv(output, sep="\t")
    assert list(counts.columns) == ["gene_id", "control_1", "treated_1"]
    assert counts.iloc[0].to_dict() == {
        "gene_id": "gene_a",
        "control_1": 7,
        "treated_1": 11,
    }


def test_nextflow_manifest_is_shiny_compatible_and_path_safe(tmp_path: Path) -> None:
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample_id\tcondition\tread1\tread2\n"
        "sample_1\tcontrol\t/private/client/R1.fastq.gz\t/private/client/R2.fastq.gz\n",
        encoding="utf-8",
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "pipeline": {"version": "1.1.0-test"},
                "project": {"id": "nextflow_test", "display_name": "Nextflow test"},
                "analysis": {"mode": "bulk_rna", "layout": "paired"},
                "rna": {"design": ["condition"]},
                "portal": {"description": "Public test fixture"},
            }
        ),
        encoding="utf-8",
    )

    fastp_dir = tmp_path / "fastp"
    fastp_dir.mkdir()
    (fastp_dir / "sample_1.fastp.json").write_text(
        json.dumps(
            {
                "summary": {
                    "before_filtering": {"total_reads": 100},
                    "after_filtering": {
                        "total_reads": 90,
                        "total_bases": 9000,
                        "q30_rate": 0.92,
                        "gc_content": 0.48,
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    (analysis_dir / "sample_1.Log.final.out").write_text(
        "Number of input reads | 90\n"
        "Uniquely mapped reads number | 80\n"
        "Uniquely mapped reads % | 88.89%\n",
        encoding="utf-8",
    )
    (analysis_dir / "gene_counts.tsv").write_text(
        "gene_id\tsample_1\ngene_a\t8\n", encoding="utf-8"
    )
    multiqc = tmp_path / "multiqc_report.html"
    multiqc.write_text("<html>test</html>\n", encoding="utf-8")
    output_root = tmp_path / "project"
    provenance_dir = tmp_path / "staged_provenance"

    build(
        argparse.Namespace(
            sample_sheet=samples,
            config=config,
            mode="bulk_rna",
            fastp_dir=fastp_dir,
            analysis_dir=analysis_dir,
            multiqc=multiqc,
            output_root=output_root,
            provenance_dir=provenance_dir,
        )
    )

    manifest_path = output_root / "portal" / "manifest.json"
    manifest_text = manifest_path.read_text(encoding="utf-8")
    manifest = json.loads(manifest_text)
    assert manifest["pipeline"]["execution_engine"] == "Nextflow"
    assert manifest["tables"]["samples"] == "portal/tables/samples.tsv"
    assert manifest["downloads"][0]["path"] == "results/multiqc/multiqc_report.html"
    assert "/private/client" not in manifest_text
    assert (provenance_dir / "config.used.yaml").is_file()
    assert (provenance_dir / "result_checksums.sha256").is_file()
