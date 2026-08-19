# Migration from the legacy draft

The uploaded legacy text was a 4,781-line pasted repository draft and ended in
the middle of an R expression-analysis module. It was not directly executable
as a repository. This release preserves its useful intent while changing the
architecture and unsafe assumptions.

| Legacy draft | Current implementation |
| --- | --- |
| `/storage/hpc/projects` | Configurable `$HOME/data`, `/cluster/pixstor`, or `/cluster/VAST` project path |
| `/scratch/hpc/temp` | `/local/scratch/$USER` through the Nextflow Hellbender profile |
| `hpc5` / `hpc5-mem` | Current public `general`; priority names stay private/configurable |
| 72-hour assembly jobs | Public jobs capped at Hellbender's current 48-hour maximum |
| Many independent array scripts | One dependency-aware, resumable Nextflow DSL2 DAG |
| One very large mutable Conda environment | Per-process environments plus a small controller environment |
| Hard-coded module versions | Current `miniconda3` bootstrap; tool versions constrained in environment files |
| BCL2FASTQ installed from Bioconda | Optional approved BCL Convert executable; no licensed redistribution |
| Generic cutadapt re-demultiplexing | Vendor dual-index demultiplexing at BCL conversion; inline barcodes are protocol-specific |
| Reference indexes built concurrently by array tasks | Single indexed dependency in the DAG or a validated reusable index |
| Hand-parsed sample rows in shell | Strict TSV/config validation before scheduling |
| Hard-coded portal usernames/passwords | Open OnDemand/Unix identity and filesystem authorization |
| Arbitrary HPC directory scanning | One configured portal root and manifest allowlist |
| Heavy DE computation inside Shiny | Batch DESeq2; Shiny only visualizes reviewed outputs |
| Incomplete/missing Shiny modules | Small complete portal with safe path and download controls |
| No automated validation | Unit tests, strict Nextflow lint, all-mode stub runs, fallback DAG dry runs, and GitHub Actions |

The release intentionally omits generic BQSR, fusion calling, ad hoc barcode
methods, and client-triggered analysis. Those require assay/reference-specific
inputs, validation, and review; presenting placeholders as an automatic full
analysis would be scientifically misleading.
