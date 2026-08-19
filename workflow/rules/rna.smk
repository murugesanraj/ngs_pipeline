"""Bulk RNA-seq alignment, quantification, tracks, and DESeq2."""


if not STATE["star_index"]:
    rule star_index:
        input:
            fasta=REF_FASTA,
            gtf=REF_GTF
        output:
            touch(STAR_INDEX_READY)
        log:
            f"{LOG_DIR}/reference/star_index.log"
        benchmark:
            f"{BENCHMARK_DIR}/reference/star_index.tsv"
        threads: 32
        resources:
            mem_mb=160000,
            runtime=720
        params:
            index_dir=STAR_INDEX_DIR,
            overhang=int(config.get("rna", {}).get("star_sjdb_overhang", 149))
        conda:
            RNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {params.index_dir:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/star_index.XXXXXX")
            cleanup() {{ rm -rf "$tmp_dir"; }}
            trap cleanup EXIT
            mkdir -p "$tmp_dir/index"
            genome_size=$(awk '!/^>/ {{n += length($0)}} END {{print n}}' {input.fasta:q})
            sa_bases=$(python -c 'import math,sys; n=max(1,int(sys.argv[1])); print(max(1,min(14,int(math.log2(n)/2-1))))' "$genome_size")
            STAR --runMode genomeGenerate \
              --runThreadN {threads} \
              --genomeDir "$tmp_dir/index" \
              --genomeFastaFiles {input.fasta:q} \
              --sjdbGTFfile {input.gtf:q} \
              --sjdbOverhang {params.overhang} \
              --genomeSAindexNbases "$sa_bases" > {log:q} 2>&1
            rm -rf {params.index_dir:q}
            mv "$tmp_dir/index" {params.index_dir:q}
            trap - EXIT
            rm -rf "$tmp_dir"
            touch {output:q}
            """


if PAIRED:
    rule star_align:
        input:
            r1=TRIM_R1_PATTERN,
            r2=TRIM_R2_PATTERN,
            index=STAR_INDEX_READY
        output:
            bam=f"{RESULTS_DIR}/rna/alignment/{{sample}}.bam",
            bai=f"{RESULTS_DIR}/rna/alignment/{{sample}}.bam.bai",
            final_log=f"{RESULTS_DIR}/rna/alignment/{{sample}}.Log.final.out",
            sj=f"{RESULTS_DIR}/rna/alignment/{{sample}}.SJ.out.tab",
            gene_counts=f"{RESULTS_DIR}/rna/alignment/{{sample}}.ReadsPerGene.out.tab"
        log:
            f"{LOG_DIR}/star/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/star/{{sample}}.tsv"
        threads: 32
        resources:
            mem_mb=64000,
            runtime=720
        params:
            index_dir=STAR_INDEX_DIR
        conda:
            RNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/star.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            prefix="$tmp_dir/{wildcards.sample}."
            STAR \
              --runThreadN {threads} \
              --genomeDir {params.index_dir:q} \
              --readFilesIn {input.r1:q} {input.r2:q} \
              --readFilesCommand zcat \
              --twopassMode Basic \
              --outSAMtype BAM SortedByCoordinate \
              --outSAMattributes NH HI AS nM MD \
              --outSAMunmapped Within \
              --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} PL:ILLUMINA \
              --quantMode GeneCounts \
              --outFileNamePrefix "$prefix" > {log:q} 2>&1
            mv "${{prefix}}Aligned.sortedByCoord.out.bam" {output.bam:q}
            mv "${{prefix}}Log.final.out" {output.final_log:q}
            mv "${{prefix}}SJ.out.tab" {output.sj:q}
            mv "${{prefix}}ReadsPerGene.out.tab" {output.gene_counts:q}
            samtools index -@ {threads} {output.bam:q} {output.bai:q}
            """
else:
    rule star_align:
        input:
            r1=TRIM_R1_PATTERN,
            index=STAR_INDEX_READY
        output:
            bam=f"{RESULTS_DIR}/rna/alignment/{{sample}}.bam",
            bai=f"{RESULTS_DIR}/rna/alignment/{{sample}}.bam.bai",
            final_log=f"{RESULTS_DIR}/rna/alignment/{{sample}}.Log.final.out",
            sj=f"{RESULTS_DIR}/rna/alignment/{{sample}}.SJ.out.tab",
            gene_counts=f"{RESULTS_DIR}/rna/alignment/{{sample}}.ReadsPerGene.out.tab"
        log:
            f"{LOG_DIR}/star/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/star/{{sample}}.tsv"
        threads: 32
        resources:
            mem_mb=64000,
            runtime=720
        params:
            index_dir=STAR_INDEX_DIR
        conda:
            RNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/star.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            prefix="$tmp_dir/{wildcards.sample}."
            STAR \
              --runThreadN {threads} \
              --genomeDir {params.index_dir:q} \
              --readFilesIn {input.r1:q} \
              --readFilesCommand zcat \
              --twopassMode Basic \
              --outSAMtype BAM SortedByCoordinate \
              --outSAMattributes NH HI AS nM MD \
              --outSAMunmapped Within \
              --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} PL:ILLUMINA \
              --quantMode GeneCounts \
              --outFileNamePrefix "$prefix" > {log:q} 2>&1
            mv "${{prefix}}Aligned.sortedByCoord.out.bam" {output.bam:q}
            mv "${{prefix}}Log.final.out" {output.final_log:q}
            mv "${{prefix}}SJ.out.tab" {output.sj:q}
            mv "${{prefix}}ReadsPerGene.out.tab" {output.gene_counts:q}
            samtools index -@ {threads} {output.bam:q} {output.bai:q}
            """


rule rna_bigwig:
    input:
        bam=f"{RESULTS_DIR}/rna/alignment/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/rna/alignment/{{sample}}.bam.bai"
    output:
        f"{RESULTS_DIR}/rna/tracks/{{sample}}.CPM.bw"
    log:
        f"{LOG_DIR}/bigwig/{{sample}}.log"
    benchmark:
        f"{BENCHMARK_DIR}/bigwig/{{sample}}.tsv"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240
    conda:
        RNA_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        bamCoverage --bam {input.bam:q} --outFileName {output:q} \
          --outFileFormat bigwig --normalizeUsing CPM --binSize 10 \
          --numberOfProcessors {threads} > {log:q} 2>&1
        """


rule featurecounts:
    input:
        bams=RNA_BAMS,
        gtf=REF_GTF
    output:
        raw=FEATURECOUNTS_RAW,
        summary=FEATURECOUNTS_SUMMARY
    log:
        f"{LOG_DIR}/featurecounts/featurecounts.log"
    benchmark:
        f"{BENCHMARK_DIR}/featurecounts/featurecounts.tsv"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=240
    params:
        strand={"unstranded": 0, "forward": 1, "reverse": 2}[
            config.get("rna", {}).get("strandedness", "unstranded")
        ],
        feature_type=config.get("rna", {}).get("featurecounts_feature_type", "exon"),
        id_attribute=config.get("rna", {}).get("featurecounts_id_attribute", "gene_id"),
        paired_flags="-p --countReadPairs -B -C" if PAIRED else ""
    conda:
        RNA_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.raw:q})" "$(dirname {log:q})"
        featureCounts -T {threads} -a {input.gtf:q} -o {output.raw:q} \
          -t {params.feature_type:q} -g {params.id_attribute:q} -s {params.strand} \
          {params.paired_flags} {input.bams:q} > {log:q} 2>&1
        test -s {output.summary:q}
        """


rule clean_featurecounts:
    input:
        FEATURECOUNTS_RAW
    output:
        COUNTS_MATRIX
    log:
        f"{LOG_DIR}/featurecounts/clean_counts.log"
    params:
        sample_ids=SAMPLE_IDS
    conda:
        RNA_ENV
    script:
        f"{SCRIPT_DIR}/aggregate_featurecounts.py"


if STATE["contrast_ids"]:
    rule deseq2:
        input:
            counts=COUNTS_MATRIX,
            samples=STATE["sample_sheet"],
            contrasts=STATE["contrast_path"]
        output:
            directory(DESEQ2_DIR)
        log:
            f"{LOG_DIR}/deseq2/deseq2.log"
        benchmark:
            f"{BENCHMARK_DIR}/deseq2/deseq2.tsv"
        resources:
            mem_mb=64000,
            runtime=360
        params:
            design=STATE["design"],
            alpha=float(config.get("rna", {}).get("alpha", 0.05)),
            min_total_count=int(config.get("rna", {}).get("min_total_count", 10))
        conda:
            REPORT_ENV
        script:
            f"{SCRIPT_DIR}/deseq2.R"
