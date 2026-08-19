# Hellbender NGS Pipeline

A reproducible Nextflow DSL2 pipeline for short-read NGS analysis on the
University of Missouri Hellbender cluster, with a read-only Shiny results portal
for authorized HPC users.

Nextflow is the primary execution engine in release 1.1.0. The repository also
keeps the validated Snakemake implementation from release 1.0.0 as a supported
fallback. Both engines use the same project YAML, sample sheet, scientific tools,
output tree, and Shiny manifest contract.

This workflow is intended for research use. It is not a validated clinical
diagnostic pipeline.

## Analysis modes

| Mode | Main outputs |
| --- | --- |
| `bulk_rna` | fastp/FastQC, STAR BAMs and junctions, featureCounts matrix, DESeq2 contrasts, PCA, sample-distance heatmap, bigWig tracks |
| `dna_reseq` | fastp/FastQC, BWA-MEM2 duplicate-marked BAMs, flagstat and mosdepth QC, optional per-sample small-variant VCFs |
| `denovo_genome` | fastp/FastQC, SPAdes or MEGAHIT contigs, QUAST, optional offline BUSCO |
| `denovo_transcriptome` | fastp/FastQC, Trinity assemblies, QUAST, optional offline BUSCO |

Every mode also creates MultiQC, checksums, Nextflow execution provenance, and
`portal/manifest.json` for the Shiny app.

## Hellbender defaults

The public configuration follows the current Hellbender documentation:

- CPU partition: `general`
- maximum public-partition wall time: 2 days
- interactive testing partition: `interactive` (maximum 4 hours)
- individual data storage: `$HOME/data`
- group storage: `/cluster/pixstor` or `/cluster/VAST`
- node-local temporary storage: `/local/scratch`
- environment module: `miniconda3`

Partition names, allocations, and module versions can change. Run the preflight
on a compute node before the first production run and check the official
[Hellbender documentation](https://docs.itrss.umsystem.edu/pub/hpc/hellbender).

## Repository layout

```text
main.nf                    Nextflow DSL2 entry point
nextflow/modules/          Common, RNA, DNA, and assembly processes
nextflow.config            Version, resources, environments, and profiles
conf/                      Hellbender, single-node, CI, and priority configs
workflow/                  Shared scripts plus the Snakemake fallback
envs/                      Pinned controller and per-process Conda environments
config/                    Safe examples only; no client data or credentials
bin/                       Bootstrap, preflight, submission, and publishing tools
slurm/                     Lightweight controller and BCL Convert jobs
app/                       Read-only Shiny results portal
docs/                      Operations, analysis, deployment, and handoff guides
tests/                     Tiny non-sensitive fixtures and automated tests
```

## Quick start on Hellbender

Do not solve environments or run analyses on a login node.

```bash
# 1. Clone into individual data or an approved group allocation.
cd "$HOME/data"
git clone https://github.com/YOUR_ORG/hellbender-ngs-pipeline.git
cd hellbender-ngs-pipeline

# 2. Obtain an interactive allocation for setup and validation.
srun -p interactive --time=02:00:00 --mem=16G \
  --cpus-per-task=4 --pty bash

# 3. Build the pinned Nextflow 26.04.6 controller and Shiny environments.
bin/bootstrap_nextflow.sh --with-shiny

# 4. Create private project files outside the tracked examples.
bin/init_project.sh "$HOME/data/ngs-projects/my_project"

# 5. Edit config.yaml, samples.tsv, and contrasts.tsv, then validate the
#    configuration and Nextflow graph while still on the compute node.
bin/preflight_nextflow.sh \
  "$HOME/data/ngs-projects/my_project/config.yaml"
```

Use absolute paths in production YAML and TSV files. The examples contain
placeholders and are intentionally rejected until copied and edited.

Submit from the login node after preflight succeeds. The wrapper itself only
calls `sbatch`; Nextflow and all analysis code run inside allocations.

```bash
bin/submit_nextflow.sh \
  --config "$HOME/data/ngs-projects/my_project/config.yaml"
```

The controller submits each process to SLURM with the CPU, memory, and time in
`nextflow.config`. It writes its resumable work directory below the project
output and records trace, report, timeline, and DAG files in `provenance/`.
Submit the same command after an interruption; `-resume` is enabled by the
wrapper and completed cached tasks are reused.

Monitor with:

```bash
squeue --me
sacct -X --starttime today \
  -o JobID,JobName%30,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

### Priority allocation

Keep allocation details out of Git:

```bash
cp conf/priority.example.config conf/priority.config
# Edit the two REPLACE_WITH_* values.

bin/submit_nextflow.sh \
  --config "$HOME/data/ngs-projects/my_project/config.yaml" \
  --nextflow-config conf/priority.config \
  --driver-partition LAB_PARTITION \
  --driver-account LAB_ACCOUNT
```

`conf/priority.config` is ignored by Git. The public `general` profile remains
the safe default.

### Single-node fallback

For a modest dataset, or if nested scheduler submission is temporarily not
appropriate at the site, request one 32-CPU/240-GB allocation and run all
processes locally inside it:

```bash
bin/submit_nextflow.sh \
  --config "$HOME/data/ngs-projects/my_project/config.yaml" \
  --single-node
```

See [docs/NEXTFLOW.md](docs/NEXTFLOW.md) for profiles, resource overrides,
resume behavior, cleanup, and direct development commands.

## Starting from BCL files

Illumina BCL Convert is licensed software and is not redistributed or installed
by this repository. If it is approved and available on Hellbender, use:

```bash
bin/submit_bcl_convert.sh \
  --run-dir /path/to/runfolder \
  --sample-sheet /path/to/SampleSheet.csv \
  --output-dir /cluster/pixstor/GROUP/project/fastq
```

Place the resulting FASTQ paths in `samples.tsv`. Standard Illumina dual indexes
should be demultiplexed by BCL Convert. Inline/UMI protocols require a dedicated,
assay-specific preprocessing module.

## Shiny results portal

The Shiny app only reads published pipeline outputs. It does not run analyses,
accept arbitrary filesystem paths, or store passwords. Each client launches it
through Hellbender Open OnDemand/RStudio, so Unix identity and approved group
permissions remain the access-control boundary.

```bash
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$HOME/data/conda/envs/hellbender-ngs-shiny"
export NGS_PORTAL_ROOT=/cluster/pixstor/GROUP/ngs_portal
Rscript app/run_app.R
```

Publish a completed project only after scientific review:

```bash
bin/publish_project.sh \
  /cluster/pixstor/GROUP/projects/my_project \
  /cluster/pixstor/GROUP/ngs_portal
```

Publishing creates a symlink and does not broaden permissions. Configure the
correct POSIX group/ACL separately under the approved data-management plan. See
[docs/PORTAL_DEPLOYMENT.md](docs/PORTAL_DEPLOYMENT.md).

## Data governance and GitHub safety

- Never commit FASTQ/BAM/VCF files, sample keys, client names, credentials,
  allocation names, or real production configuration files.
- Use coded sample identifiers and keep any re-identification key outside this
  repository and the portal.
- Cluster-attached storage is not a backup; maintain an approved archive plan.
- Confirm the permitted institutional data classification and required review
  before placing restricted data on a public or group allocation.
- Confirm University ownership and licensing requirements before publishing the
  repository under the included MIT license.

## Validation and development

With Nextflow, Java, Python, Snakemake, and ShellCheck available:

```bash
make validate
```

Focused Nextflow checks are:

```bash
nextflow lint main.nf
nextflow run . \
  -profile test \
  -params-file tests/config/bulk_rna.yaml \
  -stub-run \
  -work-dir tests/work/nextflow/bulk_rna
```

GitHub Actions performs unit tests, Python helper checks, shell lint, strict
Nextflow linting, end-to-end stub runs for every branch, and Snakemake
compatibility dry runs. Only tiny synthetic fixtures are used.

## Snakemake compatibility workflow

The release 1.0.0 interface is retained for users who need Snakemake:

```bash
bin/bootstrap_hellbender.sh --with-shiny
bin/preflight.sh --create-envs /path/to/project/config.yaml
bin/submit_pipeline.sh --config /path/to/project/config.yaml
```

Do not run both engines against the same output directory at the same time.
Choose one engine per project run; their work caches are independent.

## Documentation

- [Nextflow operation and profiles](docs/NEXTFLOW.md)
- [Hellbender setup and resource tuning](docs/HELLBENDER.md)
- [Mode-specific inputs and outputs](docs/ANALYSIS_MODES.md)
- [Secure Shiny/Open OnDemand deployment](docs/PORTAL_DEPLOYMENT.md)
- [Client handoff checklist](docs/CLIENT_HANDOFF.md)
- [Changes from the legacy scripts](docs/MIGRATION.md)

## Citation

See [CITATION.cff](CITATION.cff). Publications should also cite the primary
tools used by the selected mode and the current institutional Hellbender
acknowledgement.
