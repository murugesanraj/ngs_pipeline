from __future__ import annotations

import json
import runpy
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]


def test_manifest_contains_only_project_relative_paths(tmp_path: Path) -> None:
    project = tmp_path / "project"
    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text(
        "sample_id\tcondition\tread1\tread2\n"
        "sample_1\tcontrol\t/private/R1.fastq.gz\t/private/R2.fastq.gz\n",
        encoding="utf-8",
    )

    fastp = project / "results/qc/fastp/sample_1.fastp.json"
    fastp.parent.mkdir(parents=True)
    fastp.write_text(
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
                },
                "duplication": {"rate": 0.1},
                "adapter_cutting": {"adapter_trimmed_reads": 5},
            }
        ),
        encoding="utf-8",
    )
    star = project / "results/rna/alignment/sample_1.Log.final.out"
    star.parent.mkdir(parents=True)
    star.write_text(
        "Number of input reads | 90\n"
        "Uniquely mapped reads number | 80\n"
        "Uniquely mapped reads % | 88.89%\n",
        encoding="utf-8",
    )
    counts = project / "results/rna/counts/gene_counts.tsv"
    counts.parent.mkdir(parents=True)
    counts.write_text("gene_id\tsample_1\ngene1\t10\n", encoding="utf-8")
    multiqc = project / "results/multiqc/multiqc_report.html"
    multiqc.parent.mkdir(parents=True)
    multiqc.write_text("<html>test</html>\n", encoding="utf-8")

    outputs = SimpleNamespace(
        manifest=str(project / "portal/manifest.json"),
        samples=str(project / "portal/tables/samples.tsv"),
        qc=str(project / "portal/tables/sample_qc.tsv"),
        mapping=str(project / "portal/tables/mapping_qc.tsv"),
        config_used=str(project / "provenance/config.used.yaml"),
        checksums=str(project / "provenance/result_checksums.sha256"),
    )
    params = SimpleNamespace(
        project_dir=str(project),
        project_id="test_project",
        display_name="Test project",
        mode="bulk_rna",
        pipeline_version="test",
        sample_sheet=str(sample_sheet),
        contrast_ids=[],
        contrasts=[],
        config={"rna": {"design": ["condition"]}, "portal": {"description": "test"}},
    )
    fake = SimpleNamespace(
        params=params,
        output=outputs,
        input=SimpleNamespace(),
        log=[str(project / "logs/portal.log")],
    )
    runpy.run_path(
        str(REPO / "workflow/scripts/build_portal_manifest.py"),
        init_globals={"snakemake": fake},
    )

    manifest = json.loads(Path(outputs.manifest).read_text(encoding="utf-8"))
    assert manifest["project"]["id"] == "test_project"
    assert manifest["tables"]["samples"] == "portal/tables/samples.tsv"
    assert "/private/" not in Path(outputs.manifest).read_text(encoding="utf-8")
    assert manifest["downloads"][0]["path"] == "results/multiqc/multiqc_report.html"
