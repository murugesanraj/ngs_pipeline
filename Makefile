.PHONY: test syntax nextflow-lint nextflow-stub dry-run snakemake-dry-run validate

NEXTFLOW ?= nextflow
NEXTFLOW_CONFIGS = \
	tests/config/bulk_rna.yaml \
	tests/config/dna_reseq.yaml \
	tests/config/denovo_genome.yaml \
	tests/config/denovo_genome_megahit.yaml \
	tests/config/denovo_transcriptome.yaml

test:
	python -m pytest -q

syntax:
	python -m compileall -q workflow tests
	bash -n bin/*.sh slurm/*.sbatch

nextflow-lint:
	NXF_SYNTAX_PARSER=v2 $(NEXTFLOW) lint main.nf

nextflow-stub:
	@set -eu; for config in $(NEXTFLOW_CONFIGS); do \
		name=$$(basename "$$config" .yaml); \
		NXF_SYNTAX_PARSER=v2 $(NEXTFLOW) run . \
			-profile test -params-file "$$config" -stub-run \
			-work-dir "tests/work/nextflow/$$name"; \
	done

dry-run: nextflow-lint nextflow-stub

snakemake-dry-run:
	snakemake --snakefile workflow/Snakefile \
		--configfile tests/config/bulk_rna.yaml \
		--profile profiles/hellbender-single-node \
		--dry-run

validate: test syntax dry-run snakemake-dry-run
