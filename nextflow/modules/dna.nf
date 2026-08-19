process REFERENCE_FAIDX {
    tag 'reference'
    label 'dna'
    publishDir { "${outdir}/reference" }, mode: 'copy', overwrite: true

    input:
    path fasta
    val outdir

    output:
    path "${fasta.name}.fai", emit: index

    script:
    """
    samtools faidx ${fasta}
    """

    stub:
    """
    touch ${fasta.name}.fai
    """
}


process BWA_INDEX {
    tag 'reference'
    label 'dna'
    publishDir { "${outdir}/reference" }, mode: 'copy', overwrite: true

    input:
    path fasta
    val outdir

    output:
    path 'bwa_index', emit: index

    script:
    """
    mkdir bwa_index
    cp -L ${fasta} bwa_index/genome.fa
    bwa-mem2 index bwa_index/genome.fa
    test -s bwa_index/genome.fa.0123
    """

    stub:
    """
    mkdir bwa_index
    touch bwa_index/genome.fa
    touch bwa_index/genome.fa.0123 bwa_index/genome.fa.amb bwa_index/genome.fa.ann
    touch bwa_index/genome.fa.bwt.2bit.64 bwa_index/genome.fa.pac
    """
}


process BWA_INDEX_PREBUILT {
    tag 'prebuilt reference'
    label 'dna'

    input:
    path index_files, stageAs: 'indexes/*'
    val prefix_name
    path fasta

    output:
    path 'bwa_index', emit: index

    script:
    """
    mkdir bwa_index
    cp -L ${fasta} bwa_index/genome.fa
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
      cp -L "indexes/${prefix_name}\${suffix}" "bwa_index/genome.fa\${suffix}"
    done
    """

    stub:
    """
    mkdir bwa_index
    touch bwa_index/genome.fa
    touch bwa_index/genome.fa.0123 bwa_index/genome.fa.amb bwa_index/genome.fa.ann
    touch bwa_index/genome.fa.bwt.2bit.64 bwa_index/genome.fa.pac
    """
}


process BWA_ALIGN_MARKDUP_PE {
    tag "${meta.sample_id}"
    label 'dna'
    publishDir { "${outdir}/results/dna/alignment" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? filename : null }
    publishDir { "${outdir}/results/dna/qc" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? null : filename }

    input:
    tuple val(meta), path(read1), path(read2)
    path bwa_index
    val outdir

    output:
    tuple val(meta),
        path("${meta.sample_id}.markdup.bam"),
        path("${meta.sample_id}.markdup.bam.bai"),
        path("${meta.sample_id}.flagstat.txt"),
        path("${meta.sample_id}.idxstats.tsv"),
        path("${meta.sample_id}.markdup.stats.txt"),
        emit: alignment

    script:
    """
    bwa-mem2 mem \
      -t ${task.cpus} \
      -R '@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}\\tPL:ILLUMINA' \
      bwa_index/genome.fa ${read1} ${read2} \
      | samtools sort -n -@ ${task.cpus} -m 1G -o name.bam -
    samtools fixmate -@ ${task.cpus} -m name.bam fixmate.bam
    samtools sort -@ ${task.cpus} -m 1G -o position.bam fixmate.bam
    samtools markdup \
      -@ ${task.cpus} \
      -s \
      -f ${meta.sample_id}.markdup.stats.txt \
      position.bam ${meta.sample_id}.markdup.bam
    samtools index -@ ${task.cpus} ${meta.sample_id}.markdup.bam
    samtools flagstat -@ ${task.cpus} ${meta.sample_id}.markdup.bam > ${meta.sample_id}.flagstat.txt
    samtools idxstats -@ ${task.cpus} ${meta.sample_id}.markdup.bam > ${meta.sample_id}.idxstats.tsv
    """

    stub:
    """
    touch ${meta.sample_id}.markdup.bam ${meta.sample_id}.markdup.bam.bai
    printf '10 + 0 in total (QC-passed reads + QC-failed reads)\n9 + 0 mapped (90.00%% : N/A)\n' > ${meta.sample_id}.flagstat.txt
    touch ${meta.sample_id}.idxstats.tsv ${meta.sample_id}.markdup.stats.txt
    """
}


process BWA_ALIGN_MARKDUP_SE {
    tag "${meta.sample_id}"
    label 'dna'
    publishDir { "${outdir}/results/dna/alignment" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? filename : null }
    publishDir { "${outdir}/results/dna/qc" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? null : filename }

    input:
    tuple val(meta), path(read1)
    path bwa_index
    val outdir

    output:
    tuple val(meta),
        path("${meta.sample_id}.markdup.bam"),
        path("${meta.sample_id}.markdup.bam.bai"),
        path("${meta.sample_id}.flagstat.txt"),
        path("${meta.sample_id}.idxstats.tsv"),
        path("${meta.sample_id}.markdup.stats.txt"),
        emit: alignment

    script:
    """
    bwa-mem2 mem \
      -t ${task.cpus} \
      -R '@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}\\tPL:ILLUMINA' \
      bwa_index/genome.fa ${read1} \
      | samtools sort -n -@ ${task.cpus} -m 1G -o name.bam -
    samtools fixmate -@ ${task.cpus} -m name.bam fixmate.bam
    samtools sort -@ ${task.cpus} -m 1G -o position.bam fixmate.bam
    samtools markdup \
      -@ ${task.cpus} \
      -s \
      -f ${meta.sample_id}.markdup.stats.txt \
      position.bam ${meta.sample_id}.markdup.bam
    samtools index -@ ${task.cpus} ${meta.sample_id}.markdup.bam
    samtools flagstat -@ ${task.cpus} ${meta.sample_id}.markdup.bam > ${meta.sample_id}.flagstat.txt
    samtools idxstats -@ ${task.cpus} ${meta.sample_id}.markdup.bam > ${meta.sample_id}.idxstats.tsv
    """

    stub:
    """
    touch ${meta.sample_id}.markdup.bam ${meta.sample_id}.markdup.bam.bai
    printf '10 + 0 in total (QC-passed reads + QC-failed reads)\n9 + 0 mapped (90.00%% : N/A)\n' > ${meta.sample_id}.flagstat.txt
    touch ${meta.sample_id}.idxstats.tsv ${meta.sample_id}.markdup.stats.txt
    """
}


process MOSDEPTH {
    tag "${meta.sample_id}"
    label 'dna'
    publishDir { "${outdir}/results/dna/qc" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai)
    val outdir

    output:
    tuple val(meta),
        path("${meta.sample_id}.mosdepth.summary.txt"),
        path("${meta.sample_id}.mosdepth.global.dist.txt"),
        emit: coverage

    script:
    """
    mosdepth \
      --threads ${task.cpus} \
      --no-per-base \
      --fast-mode \
      ${meta.sample_id} ${bam}
    test -s ${meta.sample_id}.mosdepth.summary.txt
    test -s ${meta.sample_id}.mosdepth.global.dist.txt
    """

    stub:
    """
    printf 'chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t100\t100\t10.0\t0\t20\n' > ${meta.sample_id}.mosdepth.summary.txt
    printf 'total\t0\t1.0\n' > ${meta.sample_id}.mosdepth.global.dist.txt
    """
}


process CALL_VARIANTS {
    tag "${meta.sample_id}"
    label 'dna'
    publishDir { "${outdir}/results/dna/variants" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_index
    val outdir
    val minimum_quality
    val minimum_depth

    output:
    tuple val(meta),
        path("${meta.sample_id}.filtered.vcf.gz"),
        path("${meta.sample_id}.filtered.vcf.gz.tbi"),
        path("${meta.sample_id}.bcftools.stats.txt"),
        emit: results

    script:
    """
    bcftools mpileup \
      --threads ${task.cpus} \
      -f ${fasta} \
      -Ou ${bam} \
      | bcftools call \
          --threads ${task.cpus} \
          --multiallelic-caller \
          --variants-only \
          -Oz \
          -o raw.vcf.gz
    bcftools index --tbi raw.vcf.gz
    bcftools filter \
      --threads ${task.cpus} \
      --soft-filter LowQual \
      --exclude 'QUAL<${minimum_quality} || INFO/DP<${minimum_depth}' \
      -Oz \
      -o ${meta.sample_id}.filtered.vcf.gz \
      raw.vcf.gz
    bcftools index --tbi ${meta.sample_id}.filtered.vcf.gz
    bcftools stats ${meta.sample_id}.filtered.vcf.gz > ${meta.sample_id}.bcftools.stats.txt
    """

    stub:
    """
    touch ${meta.sample_id}.filtered.vcf.gz ${meta.sample_id}.filtered.vcf.gz.tbi
    printf 'SN\t0\tnumber of records:\t0\n' > ${meta.sample_id}.bcftools.stats.txt
    """
}
