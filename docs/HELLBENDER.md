# Hellbender setup and operations

This guide records the site assumptions used by release 1.1.0. Always verify
them against the live [Hellbender documentation](https://docs.itrss.umsystem.edu/pub/hpc/hellbender)
before production use.

## Storage

| Purpose | Recommended location | Notes |
| --- | --- | --- |
| User environments/small projects | `$HOME/data` | Public documentation lists 500 GB per user. Avoid the 50 GB home directory for Conda. |
| Shared production projects | `/cluster/pixstor/...` or `/cluster/VAST/...` | Use the PI's approved RDE allocation and group permissions. |
| Per-job temporary files | `/local/scratch/$USER` | Node-local, fast, and temporary. The Nextflow Hellbender profile enables process scratch here. |
| Long-term retention | Approved archive/backup | Cluster storage is not itself a backup. |

Do not hard-code a lab mount in workflow rules. Set the project output,
reference, and FASTQ paths in the private project configuration.

## Partitions and limits

The public documentation currently lists:

| Partition | Default | Maximum | Intended use |
| --- | ---: | ---: | --- |
| `general` | 1 hour | 2 days | Public CPU jobs |
| `requeue` | 10 minutes | 2 days | Preemptible/requeue public jobs |
| `interactive` | 1 hour | 4 hours | Setup, testing, and debugging |
| `gpu` | 1 hour | 2 days | Jobs that actively use GPUs |

The pipeline uses `general`; no included process requests a GPU. Check live values:

```bash
sinfo -o "%P %a %l %D %c %m"
scontrol show partition general
sinfo -o "%5D %4c %8m %28f %35G"
```

Priority partition/account names are allocation-specific. Keep them in a local,
ignored profile, not in Git.

## Login-node policy

Hellbender documentation says not to run code on login nodes. The provided
workflow wrapper uses `sbatch` to start a small controller allocation. The
controller then asks Nextflow's native SLURM executor to submit analysis jobs.
Some SLURM sites restrict nested submission; if Hellbender changes that policy,
use `--single-node` for a modest dataset or ask ITRSS for the current sanctioned
workflow-manager pattern.

Setup and preflight belong in an interactive allocation:

```bash
srun -p interactive --time=02:00:00 --mem=16G --cpus-per-task=4 --pty bash
bin/preflight_nextflow.sh /path/to/project/config.yaml
```

## Conda

The current site module name is `miniconda3`. The bootstrap stores environments
under `$HOME/data/conda/envs` by default to avoid the home quota. Nextflow process
environments are stored in `.nextflow/conda`; place the repository on data or
group storage that every compute node can read.

```bash
module load miniconda3
conda env list
```

Do not mix separately loaded Java/R/bioinformatics modules with Conda process
environments unless a tested local profile explicitly requires that setup.

## Resource tuning

Start from benchmark and SLURM accounting data rather than guesses:

```bash
sacct -X --starttime today \
  -o JobID,JobName%30,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

Override a resource without changing scientific code by creating a local,
ignored Nextflow config:

```groovy
process {
    withName: STAR_ALIGN_PE {
        memory = 96.GB
        time = 24.h
    }
}
```

Pass it with `bin/submit_nextflow.sh --nextflow-config FILE`. See
[NEXTFLOW.md](NEXTFLOW.md) for the full command and priority allocation pattern.

Keep public-partition runtime at or below 2,880 minutes. If a job reaches that
limit, first improve partitioning/checkpointing or use approved priority
resources; do not silently truncate the analysis.

## Requeue behavior

`general` is the default for reproducibility. If explicitly using `requeue`,
ensure long steps are restartable. Applications such as assemblers may not
resume inside an interrupted temporary directory. Nextflow retries common
terminated/preempted exit codes once and can reuse completed cached processes
after resubmission. Retry is not a substitute for correct memory/time requests.

## Snakemake fallback

The release 1.0.0 Snakemake profiles remain under `profiles/`. If the fallback
engine is selected, use `bin/preflight.sh`, `bin/submit_pipeline.sh`, and the
`.snakemake/conda` cache. Never run both engines against one output directory at
the same time.

## Maintenance and support

Public documentation currently schedules maintenance on the second Tuesday of
each month. Check the University status/announcement channels before long runs.
For access, allocation, network, storage, or scheduler questions, use the
current ITRSS contact information on the official Hellbender page.
