# Nextflow operation

Release 1.1.0 uses Nextflow DSL2 as its primary engine. It requires Nextflow
26.04 or later and is pinned to 26.04.6 in `envs/nextflow.yaml`. The code uses
the strict v2 syntax parser that is standard in the 26.04 release line.

## Execution model

`bin/submit_nextflow.sh` submits a small controller job. The controller runs
`slurm/nextflow_driver.sbatch`, validates the project configuration, and starts
Nextflow. Under the `hellbender` profile, each process is a separate SLURM job on
the `general` partition.

The controller must remain alive while processes are running. Hellbender's
public wall-time maximum is currently 48 hours. If a large run outlives the
controller, submit the same command again; the wrapper always uses `-resume` and
the shared work directory is preserved.

Never launch the workflow directly on a login node.

## Profiles

| Profile | Executor | Purpose |
| --- | --- | --- |
| `hellbender` | SLURM | Normal production; one scheduler job per process |
| `hellbender_single_node` | local | All processes inside one controller allocation |
| `test` | local | CI/development stub runs with resources capped at 2 CPU/4 GB |

The public profile uses `general`, process-specific CPU/memory/time requests,
`/local/scratch/$USER` for task scratch, and Conda environments cached below
`.nextflow/conda`. Keep the repository and project work directory on storage
visible from every compute node.

## Install and preflight

From an interactive compute allocation:

```bash
bin/bootstrap_nextflow.sh --with-shiny
bin/preflight_nextflow.sh /absolute/path/to/project/config.yaml
```

The bootstrap creates the controller environment under
`$HOME/data/conda/envs/hellbender-ngs-nextflow` by default. Set
`NGS_NEXTFLOW_ENV` or pass `--nextflow-prefix` to use another approved shared
location. Preflight runs the common Python validator, strict Nextflow lint, and
a graph preview without executing analysis processes.

The process environments are created lazily by Nextflow. Environment creation
happens from the allocated controller, not from the login shell, and the cache
is reused by later runs.

## Submit and resume

```bash
bin/submit_nextflow.sh --config /absolute/path/to/project/config.yaml
```

The wrapper provides:

- a stable work directory at `<output_dir>/work/nextflow`;
- automatic `-resume` behavior;
- trace, HTML report, timeline, and DAG files under `<output_dir>/provenance`;
- a controller log under `logs/nextflow-driver/` in the repository checkout;
- strict shell failure handling and private-by-default file permissions.

To resume, correct the input/configuration or resource problem and submit the
same command. Do not delete the work directory before a successful resume. A
scientific parameter change can invalidate only part of the cache, so review
the new DAG and provenance before releasing results.

Execution files include the controller job ID and SLURM restart count, for
example `nextflow_trace.123456.0.tsv`. Each resubmission therefore preserves the
earlier attempt instead of failing on or overwriting an existing HTML report.

## Priority partitions and accounts

Copy the ignored template and edit it locally:

```bash
cp conf/priority.example.config conf/priority.config
```

Then submit with both controller and task allocation details:

```bash
bin/submit_nextflow.sh \
  --config /absolute/path/to/project/config.yaml \
  --nextflow-config conf/priority.config \
  --driver-partition LAB_PARTITION \
  --driver-account LAB_ACCOUNT
```

The extra config overrides process `queue` and adds the SLURM task account. The
driver flags select the controller partition/account. Do not commit the local
file or real allocation names.

For a direct development invocation, Nextflow launcher options such as `-c`
must precede `run`:

```bash
nextflow -c conf/priority.config run . \
  -profile hellbender \
  -params-file /absolute/path/to/project/config.yaml \
  -work-dir /absolute/path/to/project/work/nextflow \
  -resume
```

Use the supplied submission wrapper for production so the controller itself is
never left on a login node.

## Resource overrides

Default resources live in `nextflow.config`. Create an ignored local config for
project/site overrides instead of editing process code:

```groovy
process {
    withName: STAR_ALIGN_PE {
        cpus = 32
        memory = 96.GB
        time = 24.h
    }
    withName: ASSEMBLE_SPADES_PE {
        memory = 220.GB
        time = 48.h
    }
}
```

Supply it with `--nextflow-config`. Keep public-partition time at or below 48
hours. Base tuning on `sacct`, the job-suffixed `nextflow_trace.*.tsv`, and the
Nextflow HTML report. More memory is not automatically better: verify node
availability and the application's own per-thread memory behavior.

## Development checks

Strict lint:

```bash
NXF_SYNTAX_PARSER=v2 nextflow lint main.nf
```

End-to-end control-flow validation without scientific tools:

```bash
nextflow run . \
  -profile test \
  -params-file tests/config/dna_reseq.yaml \
  -stub-run \
  -work-dir tests/work/nextflow/dna_reseq
```

The CI repeats stub runs for RNA, DNA, both genome assemblers, and the single-end
transcriptome path. A stub run validates scheduling and output contracts; it is
not a substitute for a small real-data acceptance run on Hellbender.

## Outputs and cleanup

Published results are copied into the configured project directory. Nextflow's
work directory contains staged inputs and intermediate products and can be much
larger. After scientific review, successful publication, and confirmation that
no resume is required, identify the exact run and preview cleanup first:

```bash
nextflow log
nextflow clean RUN_NAME -n
nextflow clean RUN_NAME -f
```

Run these commands from the same repository checkout and replace `RUN_NAME`
with the exact recorded run. Never clean a shared storage root, home directory,
or unresolved shell variable.

## Engine compatibility

Nextflow and Snakemake share project configuration and final result contracts,
but their caches are unrelated. Do not run both against one output directory at
the same time and do not expect one engine to resume the other's partial run.
Record the selected engine in the handoff; Nextflow manifests include
`pipeline.execution_engine: Nextflow`.
