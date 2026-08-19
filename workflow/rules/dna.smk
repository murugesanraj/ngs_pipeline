"""DNA resequencing alignment, coverage QC, and optional variant discovery."""


if not STATE["bwa_index_prefix"]:
    rule bwa_index:
        input:
            fasta=REF_FASTA
        output:
            touch(BWA_INDEX_READY)
        log:
            f"{LOG_DIR}/reference/bwa_index.log"
        benchmark:
            f"{BENCHMARK_DIR}/reference/bwa_index.tsv"
        threads: 16
        resources:
            mem_mb=32000,
            runtime=480
        params:
            prefix=BWA_PREFIX
        conda:
            DNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {params.prefix:q})" "$(dirname {log:q})"
            ln -sfn "$(readlink -f {input.fasta:q})" {params.prefix:q}
            bwa-mem2 index {params.prefix:q} > {log:q} 2>&1
            test -s "{params.prefix}.0123"
            touch {output:q}
            """


if PAIRED:
    rule bwa_align_markdup:
        input:
            r1=TRIM_R1_PATTERN,
            r2=TRIM_R2_PATTERN,
            index=BWA_INDEX_READY
        output:
            bam=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam",
            bai=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam.bai",
            flagstat=f"{RESULTS_DIR}/dna/qc/{{sample}}.flagstat.txt",
            idxstats=f"{RESULTS_DIR}/dna/qc/{{sample}}.idxstats.tsv",
            dup_metrics=f"{RESULTS_DIR}/dna/qc/{{sample}}.markdup.stats.txt"
        log:
            f"{LOG_DIR}/bwa/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/bwa/{{sample}}.tsv"
        threads: 32
        resources:
            mem_mb=64000,
            runtime=1440
        params:
            prefix=BWA_PREFIX
        conda:
            DNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {output.flagstat:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            exec > {log:q} 2>&1
            tmp_dir=$(mktemp -d "{resources.tmpdir}/bwa.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            rg=$(printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA' '{wildcards.sample}' '{wildcards.sample}')
            bwa-mem2 mem -t {threads} -R "$rg" {params.prefix:q} {input.r1:q} {input.r2:q} \
              | samtools sort -n -@ {threads} -m 1G -o "$tmp_dir/name.bam" -
            samtools fixmate -@ {threads} -m "$tmp_dir/name.bam" "$tmp_dir/fixmate.bam"
            samtools sort -@ {threads} -m 1G -o "$tmp_dir/position.bam" "$tmp_dir/fixmate.bam"
            samtools markdup -@ {threads} -s -f {output.dup_metrics:q} \
              "$tmp_dir/position.bam" {output.bam:q}
            samtools index -@ {threads} {output.bam:q} {output.bai:q}
            samtools flagstat -@ {threads} {output.bam:q} > {output.flagstat:q}
            samtools idxstats -@ {threads} {output.bam:q} > {output.idxstats:q}
            """
else:
    rule bwa_align_markdup:
        input:
            r1=TRIM_R1_PATTERN,
            index=BWA_INDEX_READY
        output:
            bam=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam",
            bai=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam.bai",
            flagstat=f"{RESULTS_DIR}/dna/qc/{{sample}}.flagstat.txt",
            idxstats=f"{RESULTS_DIR}/dna/qc/{{sample}}.idxstats.tsv",
            dup_metrics=f"{RESULTS_DIR}/dna/qc/{{sample}}.markdup.stats.txt"
        log:
            f"{LOG_DIR}/bwa/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/bwa/{{sample}}.tsv"
        threads: 32
        resources:
            mem_mb=64000,
            runtime=1440
        params:
            prefix=BWA_PREFIX
        conda:
            DNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {output.flagstat:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            exec > {log:q} 2>&1
            tmp_dir=$(mktemp -d "{resources.tmpdir}/bwa.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            rg=$(printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA' '{wildcards.sample}' '{wildcards.sample}')
            bwa-mem2 mem -t {threads} -R "$rg" {params.prefix:q} {input.r1:q} \
              | samtools sort -n -@ {threads} -m 1G -o "$tmp_dir/name.bam" -
            samtools fixmate -@ {threads} -m "$tmp_dir/name.bam" "$tmp_dir/fixmate.bam"
            samtools sort -@ {threads} -m 1G -o "$tmp_dir/position.bam" "$tmp_dir/fixmate.bam"
            samtools markdup -@ {threads} -s -f {output.dup_metrics:q} \
              "$tmp_dir/position.bam" {output.bam:q}
            samtools index -@ {threads} {output.bam:q} {output.bai:q}
            samtools flagstat -@ {threads} {output.bam:q} > {output.flagstat:q}
            samtools idxstats -@ {threads} {output.bam:q} > {output.idxstats:q}
            """


rule mosdepth:
    input:
        bam=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam",
        bai=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam.bai"
    output:
        summary=f"{RESULTS_DIR}/dna/qc/{{sample}}.mosdepth.summary.txt",
        distribution=f"{RESULTS_DIR}/dna/qc/{{sample}}.mosdepth.global.dist.txt"
    log:
        f"{LOG_DIR}/mosdepth/{{sample}}.log"
    benchmark:
        f"{BENCHMARK_DIR}/mosdepth/{{sample}}.tsv"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480
    params:
        prefix=lambda wc: f"{RESULTS_DIR}/dna/qc/{wc.sample}"
    conda:
        DNA_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {log:q})"
        mosdepth --threads {threads} --no-per-base --fast-mode \
          {params.prefix:q} {input.bam:q} > {log:q} 2>&1
        test -s {output.summary:q}
        test -s {output.distribution:q}
        """


if bool(config.get("dna", {}).get("call_variants", True)):
    rule call_variants:
        input:
            bam=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam",
            bai=f"{RESULTS_DIR}/dna/alignment/{{sample}}.markdup.bam.bai",
            fasta=REF_FASTA,
            fai=f"{REF_FASTA}.fai"
        output:
            vcf=f"{RESULTS_DIR}/dna/variants/{{sample}}.filtered.vcf.gz",
            tbi=f"{RESULTS_DIR}/dna/variants/{{sample}}.filtered.vcf.gz.tbi",
            stats=f"{RESULTS_DIR}/dna/variants/{{sample}}.bcftools.stats.txt"
        log:
            f"{LOG_DIR}/variants/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/variants/{{sample}}.tsv"
        threads: 8
        resources:
            mem_mb=32000,
            runtime=1440
        params:
            min_qual=float(config.get("dna", {}).get("minimum_variant_quality", 20)),
            min_depth=int(config.get("dna", {}).get("minimum_depth", 5))
        conda:
            DNA_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.vcf:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            exec > {log:q} 2>&1
            tmp_dir=$(mktemp -d "{resources.tmpdir}/variants.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            bcftools mpileup --threads {threads} -f {input.fasta:q} -Ou {input.bam:q} \
              | bcftools call --threads {threads} --multiallelic-caller --variants-only \
                  -Oz -o "$tmp_dir/raw.vcf.gz"
            bcftools index --tbi "$tmp_dir/raw.vcf.gz"
            bcftools filter --threads {threads} --soft-filter LowQual \
              --exclude 'QUAL<{params.min_qual} || INFO/DP<{params.min_depth}' \
              -Oz -o {output.vcf:q} "$tmp_dir/raw.vcf.gz"
            bcftools index --tbi -o {output.tbi:q} {output.vcf:q}
            bcftools stats {output.vcf:q} > {output.stats:q}
            """
