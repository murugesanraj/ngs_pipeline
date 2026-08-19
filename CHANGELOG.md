# Changelog

## 1.1.0 - 2026-08-19

- Added a complete Nextflow 26.04 DSL2 implementation for every analysis mode.
- Added native SLURM, single-node, test, and private priority configurations for
  Hellbender, with shared Conda caching and node-local scratch.
- Added Nextflow bootstrap, preflight, submission, resume, trace, report,
  timeline, and DAG generation.
- Kept the validated Snakemake 9 implementation as a compatibility fallback.
- Added strict Nextflow linting and end-to-end stub runs for every branch to CI.
- Added Nextflow-specific helper tests and engine provenance in Shiny manifests.

## 1.0.0 - 2026-08-19

- Replaced independent array scripts with a resumable Snakemake 9 workflow.
- Updated Hellbender defaults to `general`, `/local/scratch`, and current RDE
  mount conventions.
- Added bulk RNA-seq, DNA resequencing, de novo genome, and de novo
  transcriptome branches.
- Added per-rule Conda environments, provenance, MultiQC, tests, and CI.
- Added a read-only Shiny portal designed for Open OnDemand/Unix-group access.
- Removed hard-coded credentials and unsafe client-side analysis execution.
- Moved licensed BCL conversion to an optional, separately submitted ingress
  stage.
