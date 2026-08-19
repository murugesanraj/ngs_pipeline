process ASSEMBLE_SPADES_PE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    val outdir
    val minimum_contig_length

    output:
    tuple val(meta), path("${meta.sample_id}.contigs.fasta"), emit: contigs

    script:
    def memoryGb = task.memory.toGiga() as Integer
    """
    spades.py \
      --careful \
      --threads ${task.cpus} \
      --memory ${memoryGb} \
      -1 ${read1} \
      -2 ${read2} \
      -o spades_out
    seqkit seq -m ${minimum_contig_length} spades_out/scaffolds.fasta \
      > ${meta.sample_id}.contigs.fasta
    test -s ${meta.sample_id}.contigs.fasta
    """

    stub:
    """
    printf '>stub\nACGTACGTACGT\n' > ${meta.sample_id}.contigs.fasta
    """
}


process ASSEMBLE_SPADES_SE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1)
    val outdir
    val minimum_contig_length

    output:
    tuple val(meta), path("${meta.sample_id}.contigs.fasta"), emit: contigs

    script:
    def memoryGb = task.memory.toGiga() as Integer
    """
    spades.py \
      --careful \
      --threads ${task.cpus} \
      --memory ${memoryGb} \
      -s ${read1} \
      -o spades_out
    seqkit seq -m ${minimum_contig_length} spades_out/scaffolds.fasta \
      > ${meta.sample_id}.contigs.fasta
    test -s ${meta.sample_id}.contigs.fasta
    """

    stub:
    """
    printf '>stub\nACGTACGTACGT\n' > ${meta.sample_id}.contigs.fasta
    """
}


process ASSEMBLE_MEGAHIT_PE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    val outdir
    val minimum_contig_length

    output:
    tuple val(meta), path("${meta.sample_id}.contigs.fasta"), emit: contigs

    script:
    """
    megahit \
      -1 ${read1} \
      -2 ${read2} \
      -t ${task.cpus} \
      --min-contig-len ${minimum_contig_length} \
      -o megahit_out
    cp megahit_out/final.contigs.fa ${meta.sample_id}.contigs.fasta
    test -s ${meta.sample_id}.contigs.fasta
    """

    stub:
    """
    printf '>stub\nACGTACGTACGT\n' > ${meta.sample_id}.contigs.fasta
    """
}


process ASSEMBLE_MEGAHIT_SE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1)
    val outdir
    val minimum_contig_length

    output:
    tuple val(meta), path("${meta.sample_id}.contigs.fasta"), emit: contigs

    script:
    """
    megahit \
      -r ${read1} \
      -t ${task.cpus} \
      --min-contig-len ${minimum_contig_length} \
      -o megahit_out
    cp megahit_out/final.contigs.fa ${meta.sample_id}.contigs.fasta
    test -s ${meta.sample_id}.contigs.fasta
    """

    stub:
    """
    printf '>stub\nACGTACGTACGT\n' > ${meta.sample_id}.contigs.fasta
    """
}


process ASSEMBLE_TRINITY_PE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    val outdir
    val minimum_contig_length

    output:
    tuple val(meta), path("${meta.sample_id}.contigs.fasta"), emit: contigs

    script:
    def memoryGb = task.memory.toGiga() as Integer
    """
    Trinity \
      --seqType fq \
      --CPU ${task.cpus} \
      --max_memory ${memoryGb}G \
      --left ${read1} \
      --right ${read2} \
      --output trinity_out
    seqkit seq -m ${minimum_contig_length} trinity_out.Trinity.fasta \
      > ${meta.sample_id}.contigs.fasta
    test -s ${meta.sample_id}.contigs.fasta
    """

    stub:
    """
    printf '>stub\nACGTACGTACGT\n' > ${meta.sample_id}.contigs.fasta
    """
}


process ASSEMBLE_TRINITY_SE {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1)
    val outdir
    val minimum_contig_length

    output:
    tuple val(meta), path("${meta.sample_id}.contigs.fasta"), emit: contigs

    script:
    def memoryGb = task.memory.toGiga() as Integer
    """
    Trinity \
      --seqType fq \
      --CPU ${task.cpus} \
      --max_memory ${memoryGb}G \
      --single ${read1} \
      --output trinity_out
    seqkit seq -m ${minimum_contig_length} trinity_out.Trinity.fasta \
      > ${meta.sample_id}.contigs.fasta
    test -s ${meta.sample_id}.contigs.fasta
    """

    stub:
    """
    printf '>stub\nACGTACGTACGT\n' > ${meta.sample_id}.contigs.fasta
    """
}


process QUAST {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}/quast" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(contigs)
    val outdir
    val transcriptome_mode

    output:
    tuple val(meta),
        path("${meta.sample_id}.quast.report.tsv"),
        path("${meta.sample_id}.quast.report.html"),
        emit: reports

    script:
    def modeFlag = transcriptome_mode ? '--rna-finding' : ''
    """
    quast.py \
      --threads ${task.cpus} \
      ${modeFlag} \
      --output-dir quast_out \
      ${contigs}
    cp quast_out/report.tsv ${meta.sample_id}.quast.report.tsv
    cp quast_out/report.html ${meta.sample_id}.quast.report.html
    """

    stub:
    """
    printf 'Assembly\t# contigs\tTotal length\tN50\n${meta.sample_id}\t1\t12\t12\n' > ${meta.sample_id}.quast.report.tsv
    printf '<html><body>QUAST stub</body></html>\n' > ${meta.sample_id}.quast.report.html
    """
}


process BUSCO {
    tag "${meta.sample_id}"
    label 'assembly'
    publishDir { "${outdir}/results/assembly/${meta.sample_id}/busco" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(contigs)
    path lineage
    val outdir
    val busco_mode

    output:
    tuple val(meta), path("${meta.sample_id}.busco"), emit: results

    script:
    """
    busco \
      --offline \
      --cpu ${task.cpus} \
      --mode ${busco_mode} \
      --lineage_dataset ${lineage} \
      --in ${contigs} \
      --out ${meta.sample_id}.busco
    """

    stub:
    """
    mkdir ${meta.sample_id}.busco
    touch ${meta.sample_id}.busco/short_summary.stub.txt
    """
}
