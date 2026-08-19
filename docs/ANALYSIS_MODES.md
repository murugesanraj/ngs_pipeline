# Analysis modes

All paths in production configuration should be absolute. All samples in one
run must use the same read layout. Sample and contrast identifiers may contain
letters, digits, dots, underscores, and hyphens only.

## Common input

`samples.tsv` requires:

| Column | Requirement |
| --- | --- |
| `sample_id` | Required, unique, coded identifier |
| `read1` | Required readable FASTQ/FASTQ.GZ path |
| `read2` | Required for `paired`; blank/absent for `single` |
| RNA design columns | Each term listed in `rna.design` must exist and be non-empty |

FASTQ integrity should be verified against instrument/vendor checksums before
starting. The pipeline preserves raw inputs and writes trimmed reads under the
project results directory.

## `bulk_rna`

Required reference inputs:

- genome FASTA;
- matching gene annotation GTF;
- optional prebuilt STAR index made from the same FASTA/GTF/read-length model.

The workflow performs STAR two-pass alignment, featureCounts exon-to-gene
aggregation, CPM-normalized bigWig generation, and DESeq2. `rna.design` is an
ordered list of categorical covariates. Contrasts are explicit:

```text
contrast_id  factor     numerator  denominator  label
drug_vs_ctl  condition  drug       control      Drug vs control
```

The numerator/denominator direction controls the log2 fold-change sign. The app
does not recompute DE results. Review count filtering, model rank, sample
replication, outliers, hidden confounding, and reference/annotation compatibility
before interpreting a contrast.

Important limitations:

- gene-level counts only; no transcript-level inference;
- no automatic surrogate-variable or unwanted-variation correction;
- no automatic gene-symbol annotation or enrichment analysis;
- factor terms are treated as categorical;
- fusion detection from the legacy draft is intentionally omitted because a
  trustworthy fusion workflow requires dedicated callers, filtering, and
  validation.

## `dna_reseq`

Required reference input: genome FASTA, with an optional matching BWA-MEM2 index.

The branch uses BWA-MEM2, name/coordinate sorting, `samtools fixmate` and
`markdup`, flagstat, idxstats, and mosdepth. Optional per-sample variants use
`bcftools mpileup/call` and a clearly labeled soft filter.

Important limitations:

- the included VCF step is per-sample research discovery, not joint cohort
  genotyping;
- no base-quality recalibration, VQSR, pedigree model, CNV/SV caller, somatic
  model, or clinical interpretation;
- targeted panels/exomes need capture intervals and coverage thresholds not
  supplied by this generic mode;
- use a validated, assay-specific workflow for regulated or diagnostic work.

Set `dna.call_variants: false` when only analysis-ready BAM/coverage QC is
required.

## `denovo_genome`

No reference is required. Select SPAdes for isolate/small-genome projects or
MEGAHIT for large/metagenomic read sets. Assemblies are per sample, followed by
QUAST. Optional BUSCO uses a local lineage directory with `--offline`.

Short-read assembly alone may not resolve repeats, heterozygosity, haplotypes,
or chromosome structure. Carefully review genome size, ploidy, contamination,
coverage, and the biological appropriateness of per-sample assembly.

## `denovo_transcriptome`

No reference is required. Trinity assembles each sample independently, followed
by QUAST and optional offline BUSCO. A biological study may instead need pooled
assembly, normalization, redundancy reduction, coding-sequence prediction,
contaminant screening, annotation, and abundance-aware filtering. Those choices
are organism- and study-specific and are not silently applied here.

## Output contract

Each project contains:

```text
results/       Scientific outputs and MultiQC
reference/     Project-local reference links/indexes when needed
provenance/    Config, checksums, Nextflow trace/report/timeline/DAG
portal/        Sanitized tables and manifest consumed by Shiny
work/nextflow/ Resumable task cache and staged intermediates
```

The private `provenance/config.used.yaml` contains source paths and should not be
published outside the authorized project team. The portal manifest stores only
project-relative result paths. The controller log is written below
`logs/nextflow-driver/` in the repository checkout. The Snakemake fallback keeps
its rule logs and benchmarks under the configured project output.
