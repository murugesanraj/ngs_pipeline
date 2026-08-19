from __future__ import annotations

import copy
import csv
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflow" / "lib"))

from pipeline_config import ConfigError, validate_and_load  # noqa: E402


@pytest.fixture()
def config() -> dict:
    with (REPO / "tests" / "config" / "bulk_rna.yaml").open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    loaded["project"]["output_dir"] = str(REPO / "tests" / "work" / "bulk_rna")
    loaded["project"]["sample_sheet"] = str(REPO / "tests" / "config" / "samples.tsv")
    loaded["reference"]["fasta"] = str(REPO / "tests" / "data" / "genome.fa")
    loaded["reference"]["annotation_gtf"] = str(REPO / "tests" / "data" / "genes.gtf")
    loaded["rna"]["contrasts"] = str(REPO / "tests" / "config" / "contrasts.tsv")
    return loaded


def test_valid_bulk_rna_config(config: dict, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(REPO)
    state = validate_and_load(config)
    assert state["mode"] == "bulk_rna"
    assert state["paired"] is True
    assert state["sample_ids"] == ["control_1", "control_2", "treated_1", "treated_2"]
    assert state["contrast_ids"] == ["treated_vs_control"]


def test_rejects_unsafe_project_identifier(config: dict) -> None:
    bad = copy.deepcopy(config)
    bad["project"]["id"] = "../another_project"
    with pytest.raises(ConfigError, match="project.id"):
        validate_and_load(bad)


def test_rejects_placeholder_output(config: dict) -> None:
    bad = copy.deepcopy(config)
    bad["project"]["output_dir"] = "/cluster/pixstor/GROUP/project"
    with pytest.raises(ConfigError, match="placeholder"):
        validate_and_load(bad)


def test_rejects_missing_paired_read(config: dict, tmp_path: Path) -> None:
    source = REPO / "tests" / "config" / "samples.tsv"
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[0]["read2"] = ""
    sheet = tmp_path / "samples.tsv"
    with sheet.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    bad = copy.deepcopy(config)
    bad["project"]["sample_sheet"] = str(sheet)
    with pytest.raises(ConfigError, match="has no read2"):
        validate_and_load(bad)


def test_rejects_unknown_contrast_group(config: dict, tmp_path: Path) -> None:
    contrast = tmp_path / "contrasts.tsv"
    contrast.write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\tlabel\n"
        "bad\tcondition\tmissing\tcontrol\tBad\n",
        encoding="utf-8",
    )
    bad = copy.deepcopy(config)
    bad["rna"]["contrasts"] = str(contrast)
    with pytest.raises(ConfigError, match="numerator"):
        validate_and_load(bad)
