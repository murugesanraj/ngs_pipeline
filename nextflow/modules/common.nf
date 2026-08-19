process FASTQC_RAW {
    tag "${meta.sample_id}:${read_label}"
    label 'qc'
    publishDir { "${outdir}/results/qc/fastqc/raw" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), val(read_label), path(read_file)
    val outdir

    output:
    tuple val(meta), val(read_label),
        path("${meta.sample_id}.${read_label}.raw_fastqc.html"),
        path("${meta.sample_id}.${read_label}.raw_fastqc.zip"),
        emit: reports

    script:
    """
    mkdir fastqc_tmp
    fastqc --threads ${task.cpus} --outdir fastqc_tmp ${read_file}
    mv fastqc_tmp/*_fastqc.html ${meta.sample_id}.${read_label}.raw_fastqc.html
    mv fastqc_tmp/*_fastqc.zip ${meta.sample_id}.${read_label}.raw_fastqc.zip
    """

    stub:
    """
    touch ${meta.sample_id}.${read_label}.raw_fastqc.html
    touch ${meta.sample_id}.${read_label}.raw_fastqc.zip
    """
}


process FASTP_PE {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir { "${outdir}/results/qc/fastp" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    val outdir
    val qualified_quality
    val unqualified_limit
    val minimum_length
    val detect_adapters

    output:
    tuple val(meta),
        path("${meta.sample_id}.R1.trimmed.fastq.gz"),
        path("${meta.sample_id}.R2.trimmed.fastq.gz"),
        emit: reads
    tuple val(meta),
        path("${meta.sample_id}.fastp.json"),
        path("${meta.sample_id}.fastp.html"),
        emit: reports

    script:
    def detectFlag = detect_adapters ? '--detect_adapter_for_pe' : ''
    """
    fastp \
      --in1 ${read1} --in2 ${read2} \
      --out1 ${meta.sample_id}.R1.trimmed.fastq.gz \
      --out2 ${meta.sample_id}.R2.trimmed.fastq.gz \
      --json ${meta.sample_id}.fastp.json \
      --html ${meta.sample_id}.fastp.html \
      --thread ${task.cpus} \
      --qualified_quality_phred ${qualified_quality} \
      --unqualified_percent_limit ${unqualified_limit} \
      --length_required ${minimum_length} \
      --correction ${detectFlag}
    """

    stub:
    """
    touch ${meta.sample_id}.R1.trimmed.fastq.gz
    touch ${meta.sample_id}.R2.trimmed.fastq.gz
    printf '{"summary":{"before_filtering":{"total_reads":10},"after_filtering":{"total_reads":9,"total_bases":900,"q30_rate":0.9,"gc_content":0.5}}}\n' > ${meta.sample_id}.fastp.json
    touch ${meta.sample_id}.fastp.html
    """
}


process FASTP_SE {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir { "${outdir}/results/qc/fastp" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1)
    val outdir
    val qualified_quality
    val unqualified_limit
    val minimum_length

    output:
    tuple val(meta), path("${meta.sample_id}.R1.trimmed.fastq.gz"), emit: reads
    tuple val(meta),
        path("${meta.sample_id}.fastp.json"),
        path("${meta.sample_id}.fastp.html"),
        emit: reports

    script:
    """
    fastp \
      --in1 ${read1} \
      --out1 ${meta.sample_id}.R1.trimmed.fastq.gz \
      --json ${meta.sample_id}.fastp.json \
      --html ${meta.sample_id}.fastp.html \
      --thread ${task.cpus} \
      --qualified_quality_phred ${qualified_quality} \
      --unqualified_percent_limit ${unqualified_limit} \
      --length_required ${minimum_length}
    """

    stub:
    """
    touch ${meta.sample_id}.R1.trimmed.fastq.gz
    printf '{"summary":{"before_filtering":{"total_reads":10},"after_filtering":{"total_reads":9,"total_bases":900,"q30_rate":0.9,"gc_content":0.5}}}\n' > ${meta.sample_id}.fastp.json
    touch ${meta.sample_id}.fastp.html
    """
}


process FASTQC_TRIMMED {
    tag "${meta.sample_id}:${read_label}"
    label 'qc'
    publishDir { "${outdir}/results/qc/fastqc/trimmed" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), val(read_label), path(read_file)
    val outdir

    output:
    tuple val(meta), val(read_label),
        path("${meta.sample_id}.${read_label}.trimmed_fastqc.html"),
        path("${meta.sample_id}.${read_label}.trimmed_fastqc.zip"),
        emit: reports

    script:
    """
    mkdir fastqc_tmp
    fastqc --threads ${task.cpus} --outdir fastqc_tmp ${read_file}
    mv fastqc_tmp/*_fastqc.html ${meta.sample_id}.${read_label}.trimmed_fastqc.html
    mv fastqc_tmp/*_fastqc.zip ${meta.sample_id}.${read_label}.trimmed_fastqc.zip
    """

    stub:
    """
    touch ${meta.sample_id}.${read_label}.trimmed_fastqc.html
    touch ${meta.sample_id}.${read_label}.trimmed_fastqc.zip
    """
}


process FASTQ_SCREEN {
    tag "${meta.sample_id}"
    label 'qc'
    publishDir { "${outdir}/results/qc/fastq_screen" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1)
    path screen_config
    val outdir
    val subset_size

    output:
    tuple val(meta),
        path("${meta.sample_id}_screen.txt"),
        path("${meta.sample_id}_screen.html"),
        emit: reports

    script:
    """
    mkdir screen_tmp
    fastq_screen \
      --conf ${screen_config} \
      --subset ${subset_size} \
      --threads ${task.cpus} \
      --outdir screen_tmp \
      ${read1}
    mv screen_tmp/*_screen.txt ${meta.sample_id}_screen.txt
    mv screen_tmp/*_screen.html ${meta.sample_id}_screen.html
    """

    stub:
    """
    touch ${meta.sample_id}_screen.txt
    touch ${meta.sample_id}_screen.html
    """
}


process MULTIQC {
    tag 'all samples'
    label 'qc'
    publishDir { "${outdir}/results/multiqc" }, mode: 'copy', overwrite: true

    input:
    path qc_files, stageAs: 'inputs/*'
    val outdir

    output:
    path 'multiqc_report.html', emit: report

    script:
    """
    multiqc --force --filename multiqc_report.html inputs
    """

    stub:
    """
    printf '<html><body>Nextflow stub MultiQC report</body></html>\n' > multiqc_report.html
    """
}


process BUILD_PORTAL {
    tag 'client manifest'
    label 'report'
    publishDir { outdir }, mode: 'copy', overwrite: true, pattern: 'portal'
    publishDir { "${outdir}/provenance" }, mode: 'copy', overwrite: true,
        pattern: 'config.used.yaml'
    publishDir { "${outdir}/provenance" }, mode: 'copy', overwrite: true,
        pattern: 'result_checksums.sha256'

    input:
    path sample_sheet, stageAs: 'metadata/samples.tsv'
    path project_config, stageAs: 'metadata/config.json'
    val analysis_mode
    path fastp_jsons, stageAs: 'fastp/*'
    path analysis_files, stageAs: 'analysis/*'
    path multiqc_report, stageAs: 'multiqc/multiqc_report.html'
    path portal_script, stageAs: 'scripts/build_portal_manifest_nextflow.py'
    val outdir

    output:
    path 'portal', emit: portal
    path 'config.used.yaml', emit: config_used
    path 'result_checksums.sha256', emit: checksums

    script:
    """
    python scripts/build_portal_manifest_nextflow.py \
      --sample-sheet metadata/samples.tsv \
      --config metadata/config.json \
      --mode ${analysis_mode} \
      --fastp-dir fastp \
      --analysis-dir analysis \
      --multiqc multiqc/multiqc_report.html \
      --output-root . \
      --provenance-dir .
    """

    stub:
    """
    mkdir -p portal/tables
    printf 'sample_id\n' > portal/tables/samples.tsv
    printf 'sample_id\n' > portal/tables/sample_qc.tsv
    printf 'sample_id\n' > portal/tables/mapping_qc.tsv
    printf '{"schema_version":"1.0","project":{"id":"stub"}}\n' > portal/manifest.json
    printf 'stub: true\n' > config.used.yaml
    touch result_checksums.sha256
    """
}
