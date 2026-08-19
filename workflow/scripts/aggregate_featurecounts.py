"""Convert featureCounts output into a clean integer gene-by-sample matrix."""

from pathlib import Path

import pandas as pd

raw_path = Path(str(snakemake.input[0]))  # noqa: F821
output_path = Path(str(snakemake.output[0]))  # noqa: F821
log_path = Path(str(snakemake.log[0]))  # noqa: F821
sample_ids = list(snakemake.params.sample_ids)  # noqa: F821

output_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

frame = pd.read_csv(raw_path, sep="\t", comment="#")
annotation_columns = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
missing = [column for column in annotation_columns if column not in frame.columns]
if missing:
    raise ValueError(f"featureCounts output is missing columns: {', '.join(missing)}")

count_columns = [column for column in frame.columns if column not in annotation_columns]
if len(count_columns) != len(sample_ids):
    raise ValueError(
        f"featureCounts produced {len(count_columns)} count columns for {len(sample_ids)} samples"
    )

clean = frame[["Geneid", *count_columns]].copy()
clean.columns = ["gene_id", *sample_ids]
for sample_id in sample_ids:
    clean[sample_id] = pd.to_numeric(clean[sample_id], errors="raise").astype("int64")

if clean["gene_id"].duplicated().any():
    duplicates = clean.loc[clean["gene_id"].duplicated(), "gene_id"].head(5).tolist()
    raise ValueError(f"Duplicate gene identifiers in featureCounts output: {duplicates}")

clean.to_csv(output_path, sep="\t", index=False)
log_path.write_text(
    f"genes={len(clean)}\nsamples={len(sample_ids)}\nsource={raw_path}\n",
    encoding="utf-8",
)
