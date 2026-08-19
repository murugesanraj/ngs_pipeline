# Contributing

1. Create a branch from the default branch.
2. Do not use production or identifiable client data in tests or examples.
3. Keep cluster-specific settings in profiles and scientific settings in the
   project configuration.
4. Add or update tests for configuration, path handling, and DAG construction.
5. Run `make validate` before opening a pull request. At minimum, run
   `make test`, `make nextflow-lint`, and the stub workflow for the mode changed.
6. Describe scientific behavior changes and any migration requirement.

Tool upgrades should be made deliberately in `envs/*.yaml`, tested on a small
reference dataset, and recorded in `CHANGELOG.md`. Do not silently change a
reference build, genome annotation, counting model, or statistical design.

Nextflow is the primary engine. Keep `main.nf`, `nextflow/modules/`, and the
Hellbender configs in `conf/` compatible with strict DSL2 syntax. The Snakemake
implementation remains a supported fallback, so shared scientific behavior and
the Shiny manifest contract must stay aligned across both engines.
