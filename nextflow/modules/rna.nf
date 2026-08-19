process STAR_INDEX {
    tag 'reference'
    label 'rna'
    publishDir { "${outdir}/reference" }, mode: 'copy', overwrite: true

    input:
    path fasta
    path gtf
    val outdir
    val sjdb_overhang

    output:
    path 'star_index', emit: index

    script:
    """
    mkdir star_index
    genome_size=\$(awk '!/^>/ {n += length(\$0)} END {print n}' ${fasta})
    sa_bases=\$(python -c 'import math,sys; n=max(1,int(sys.argv[1])); print(max(1,min(14,int(math.log2(n)/2-1))))' "\$genome_size")
    STAR \
      --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir star_index \
      --genomeFastaFiles ${fasta} \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang ${sjdb_overhang} \
      --genomeSAindexNbases "\$sa_bases"
    touch star_index/.complete
    """

    stub:
    """
    mkdir star_index
    touch star_index/Genome star_index/SA star_index/SAindex
    touch star_index/chrLength.txt star_index/chrName.txt star_index/.complete
    """
}


process STAR_ALIGN_PE {
    tag "${meta.sample_id}"
    label 'rna'
    publishDir { "${outdir}/results/rna/alignment" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    path star_index
    val outdir

    output:
    tuple val(meta),
        path("${meta.sample_id}.bam"),
        path("${meta.sample_id}.bam.bai"),
        path("${meta.sample_id}.Log.final.out"),
        path("${meta.sample_id}.SJ.out.tab"),
        path("${meta.sample_id}.ReadsPerGene.out.tab"),
        emit: alignment

    script:
    """
    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${star_index} \
      --readFilesIn ${read1} ${read2} \
      --readFilesCommand zcat \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --outSAMunmapped Within \
      --outSAMattrRGline ID:${meta.sample_id} SM:${meta.sample_id} PL:ILLUMINA \
      --quantMode GeneCounts \
      --outFileNamePrefix ${meta.sample_id}.
    mv ${meta.sample_id}.Aligned.sortedByCoord.out.bam ${meta.sample_id}.bam
    samtools index -@ ${task.cpus} ${meta.sample_id}.bam
    """

    stub:
    """
    touch ${meta.sample_id}.bam ${meta.sample_id}.bam.bai
    printf 'Number of input reads | 10\nUniquely mapped reads number | 9\nUniquely mapped reads %% | 90.0%%\n' > ${meta.sample_id}.Log.final.out
    touch ${meta.sample_id}.SJ.out.tab ${meta.sample_id}.ReadsPerGene.out.tab
    """
}


process STAR_ALIGN_SE {
    tag "${meta.sample_id}"
    label 'rna'
    publishDir { "${outdir}/results/rna/alignment" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1)
    path star_index
    val outdir

    output:
    tuple val(meta),
        path("${meta.sample_id}.bam"),
        path("${meta.sample_id}.bam.bai"),
        path("${meta.sample_id}.Log.final.out"),
        path("${meta.sample_id}.SJ.out.tab"),
        path("${meta.sample_id}.ReadsPerGene.out.tab"),
        emit: alignment

    script:
    """
    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir ${star_index} \
      --readFilesIn ${read1} \
      --readFilesCommand zcat \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --outSAMunmapped Within \
      --outSAMattrRGline ID:${meta.sample_id} SM:${meta.sample_id} PL:ILLUMINA \
      --quantMode GeneCounts \
      --outFileNamePrefix ${meta.sample_id}.
    mv ${meta.sample_id}.Aligned.sortedByCoord.out.bam ${meta.sample_id}.bam
    samtools index -@ ${task.cpus} ${meta.sample_id}.bam
    """

    stub:
    """
    touch ${meta.sample_id}.bam ${meta.sample_id}.bam.bai
    printf 'Number of input reads | 10\nUniquely mapped reads number | 9\nUniquely mapped reads %% | 90.0%%\n' > ${meta.sample_id}.Log.final.out
    touch ${meta.sample_id}.SJ.out.tab ${meta.sample_id}.ReadsPerGene.out.tab
    """
}


process RNA_BIGWIG {
    tag "${meta.sample_id}"
    label 'rna'
    publishDir { "${outdir}/results/rna/tracks" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai)
    val outdir

    output:
    tuple val(meta), path("${meta.sample_id}.CPM.bw"), emit: tracks

    script:
    """
    bamCoverage \
      --bam ${bam} \
      --outFileName ${meta.sample_id}.CPM.bw \
      --outFileFormat bigwig \
      --normalizeUsing CPM \
      --binSize 10 \
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch ${meta.sample_id}.CPM.bw
    """
}


process FEATURECOUNTS {
    tag 'all samples'
    label 'rna'
    publishDir { "${outdir}/results/rna/counts" }, mode: 'copy', overwrite: true

    input:
    path bams, stageAs: 'bams/*'
    path gtf
    val outdir
    val paired
    val strandedness
    val feature_type
    val id_attribute

    output:
    path 'featurecounts.raw.tsv', emit: raw
    path 'featurecounts.raw.tsv.summary', emit: summary

    script:
    def strandCode = ['unstranded': 0, 'forward': 1, 'reverse': 2][strandedness]
    if (strandCode == null) {
        error "Unsupported RNA strandedness: ${strandedness}"
    }
    def pairedFlags = paired ? '-p --countReadPairs -B -C' : ''
    """
    featureCounts \
      -T ${task.cpus} \
      -a ${gtf} \
      -o featurecounts.raw.tsv \
      -t ${feature_type} \
      -g ${id_attribute} \
      -s ${strandCode} \
      ${pairedFlags} \
      bams/*.bam
    test -s featurecounts.raw.tsv.summary
    """

    stub:
    """
    printf 'Geneid\tChr\tStart\tEnd\tStrand\tLength\tstub.bam\ngene1\tchr1\t1\t10\t+\t10\t5\n' > featurecounts.raw.tsv
    printf 'Status\tstub.bam\nAssigned\t5\n' > featurecounts.raw.tsv.summary
    """
}


process CLEAN_FEATURECOUNTS {
    tag 'count matrix'
    label 'rna'
    publishDir { "${outdir}/results/rna/counts" }, mode: 'copy', overwrite: true

    input:
    path raw_counts
    path clean_script, stageAs: 'aggregate_featurecounts_nextflow.py'
    val outdir

    output:
    path 'gene_counts.tsv', emit: counts

    script:
    """
    python aggregate_featurecounts_nextflow.py \
      --input ${raw_counts} \
      --output gene_counts.tsv
    """

    stub:
    """
    printf 'gene_id\tstub\ngene1\t5\n' > gene_counts.tsv
    """
}


process DESEQ2 {
    tag 'all contrasts'
    label 'report'
    publishDir { "${outdir}/results/rna" }, mode: 'copy', overwrite: true

    input:
    path counts
    path sample_sheet, stageAs: 'samples.tsv'
    path contrasts, stageAs: 'contrasts.tsv'
    path deseq2_script, stageAs: 'deseq2_nextflow.R'
    val outdir
    env 'DESEQ2_DESIGN'
    val alpha
    val minimum_total_count

    output:
    path 'deseq2', emit: results

    script:
    """
    Rscript deseq2_nextflow.R \
      --counts ${counts} \
      --samples samples.tsv \
      --contrasts contrasts.tsv \
      --output deseq2 \
      --alpha ${alpha} \
      --min-total-count ${minimum_total_count}
    """

    stub:
    """
    mkdir -p deseq2/stub_contrast
    printf 'gene_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\ngene1\t5\t1\t0.1\t0.2\n' > deseq2/stub_contrast/results.tsv
    printf 'contrast_id\tlabel\talpha\nstub_contrast\tStub contrast\t0.05\n' > deseq2/contrast_summary.tsv
    printf 'gene_id\tstub\ngene1\t5\n' > deseq2/normalized_counts.tsv
    touch deseq2/pca.png deseq2/sample_distance.png
    touch deseq2/stub_contrast/volcano.png deseq2/stub_contrast/ma.png
    """
}
