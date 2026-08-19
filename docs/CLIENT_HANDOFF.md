# Client handoff checklist

## Before analysis

- Confirm scope, organism/build, annotation release, assay/library kit, read
  layout, strandedness, expected comparisons, and deliverables.
- Confirm coded sample identifiers and a separate owner for any re-identification
  key.
- Confirm RDE location, access group, retention, backup/archive, and data class.
- Verify vendor checksums and sample-sheet/index assignments.
- Freeze the project configuration in version control without client data.

## Before release

- Review FastQC/fastp and contamination screening when configured.
- Review mapping/coverage or assembly completeness in context; generic thresholds
  are not automatic biological pass/fail criteria.
- For RNA-seq, inspect PCA/sample distance, design rank, replicate structure,
  contrast direction, count filtering, and influential samples.
- For DNA, state whether VCFs are exploratory per-sample calls and list omitted
  clinical/cohort processing.
- Review all portal-facing labels and notes for privacy.
- Confirm reference files and software environment records.
- Record exceptions, reruns, and manual decisions.

## Handoff package

- MultiQC HTML report;
- sanitized sample/QC tables;
- mode-specific result tables and figures;
- pipeline release/commit, execution engine, resolved config, and
  reference/annotation versions;
- interpretation/limitations note;
- retention and support contact;
- checksum file for client-facing small outputs.

Do not hand off only the Shiny URL. The portal is a convenience view; the
versioned result package and analysis record remain the reproducible deliverable.
