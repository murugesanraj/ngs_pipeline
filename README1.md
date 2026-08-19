# ngs_pipeline
End-to-end parallelized Next-Generation Sequencing (NGS) analysis pipeline optimized for the HPC cluster, running from raw BCL conversion to multi-omics assembly, mapping, and RStudio reporting.

## Pipeline Workflow
1. **Base Calling & Demultiplexing:** `bcl2fastq` conversion and barcode splitting.
2. **Quality Control & Contamination Screening:** `fastp` for adapter/quality trimming and `fastq_screen` for contaminant checks.
3. **Alignment & Assembly:** 
   - *De Novo* Assembly (Genome / Transcriptome)
   - Reference-based Genome Mapping
   - Reference-based Transcriptome Mapping using `STAR`
4. **Downstream Reporting:** Processed count matrices and summaries structured for RStudio Server visualization.
