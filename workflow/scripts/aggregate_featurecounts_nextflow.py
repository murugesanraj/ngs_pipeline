#!/usr/bin/env python3
"""Convert Nextflow featureCounts output to a clean integer count matrix."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,99}$")
ANNOTATION_COLUMNS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]


def sample_id_from_column(column: str) -> str:
    name = Path(column).name
    if name.endswith(".markdup.bam"):
        name = name.removesuffix(".markdup.bam")
    elif name.endswith(".bam"):
        name = name.removesuffix(".bam")
    if not SAFE_ID.fullmatch(name):
        raise ValueError(f"Unsafe sample identifier derived from featureCounts: {name!r}")
    return name


def clean_featurecounts(input_path: Path, output_path: Path) -> None:
    frame = pd.read_csv(input_path, sep="\t", comment="#")
    missing = [column for column in ANNOTATION_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"featureCounts output is missing columns: {', '.join(missing)}")

    count_columns = [column for column in frame.columns if column not in ANNOTATION_COLUMNS]
    if not count_columns:
        raise ValueError("featureCounts output has no sample columns")
    sample_ids = [sample_id_from_column(column) for column in count_columns]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("featureCounts column names produce duplicate sample identifiers")

    clean = frame[["Geneid", *count_columns]].copy()
    clean.columns = ["gene_id", *sample_ids]
    for sample_id in sample_ids:
        clean[sample_id] = pd.to_numeric(clean[sample_id], errors="raise").astype("int64")
    if clean["gene_id"].duplicated().any():
        duplicates = clean.loc[clean["gene_id"].duplicated(), "gene_id"].head(5).tolist()
        raise ValueError(f"Duplicate gene identifiers: {duplicates}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    clean.to_csv(output_path, sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    clean_featurecounts(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
