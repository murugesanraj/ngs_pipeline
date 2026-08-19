"""Common FASTQ QC, reference staging, MultiQC, and portal rules."""


def raw_read(wildcards):
    key = "read1" if wildcards.read == "R1" else "read2"
    value = SAMPLES[wildcards.sample][key]
    if not value:
        raise ValueError(f"No {key} configured for {wildcards.sample}")
    return value


def trimmed_read(wildcards):
    if wildcards.read == "R1":
        return TRIM_R1_PATTERN.format(sample=wildcards.sample)
    return TRIM_R2_PATTERN.format(sample=wildcards.sample)


rule fastqc_raw:
    input:
        raw_read
    output:
        html=f"{RESULTS_DIR}/qc/fastqc/raw/{{sample}}.{{read}}_fastqc.html",
        zip=f"{RESULTS_DIR}/qc/fastqc/raw/{{sample}}.{{read}}_fastqc.zip"
    log:
        f"{LOG_DIR}/fastqc/raw/{{sample}}.{{read}}.log"
    benchmark:
        f"{BENCHMARK_DIR}/fastqc/raw/{{sample}}.{{read}}.tsv"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        QC_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.html:q})" "$(dirname {log:q})" "{resources.tmpdir}"
        tmp_dir=$(mktemp -d "{resources.tmpdir}/fastqc.raw.XXXXXX")
        trap 'rm -rf "$tmp_dir"' EXIT
        fastqc --threads {threads} --outdir "$tmp_dir" {input:q} > {log:q} 2>&1
        mv "$tmp_dir"/*_fastqc.html {output.html:q}
        mv "$tmp_dir"/*_fastqc.zip {output.zip:q}
        """


if PAIRED:
    rule fastp:
        input:
            r1=lambda wc: SAMPLES[wc.sample]["read1"],
            r2=lambda wc: SAMPLES[wc.sample]["read2"]
        output:
            r1=TRIM_R1_PATTERN,
            r2=TRIM_R2_PATTERN,
            json=FASTP_JSON_PATTERN,
            html=FASTP_HTML_PATTERN
        log:
            f"{LOG_DIR}/fastp/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/fastp/{{sample}}.tsv"
        threads: 16
        resources:
            mem_mb=16000,
            runtime=240
        params:
            q=config.get("qc", {}).get("fastp", {}).get("qualified_quality_phred", 20),
            unqualified=config.get("qc", {}).get("fastp", {}).get("unqualified_percent_limit", 40),
            min_length=config.get("qc", {}).get("fastp", {}).get("length_required", 36),
            detect=(
                "--detect_adapter_for_pe"
                if config.get("qc", {}).get("fastp", {}).get("detect_adapter_for_pe", True)
                else ""
            )
        conda:
            QC_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"
            fastp \
              --in1 {input.r1:q} --in2 {input.r2:q} \
              --out1 {output.r1:q} --out2 {output.r2:q} \
              --json {output.json:q} --html {output.html:q} \
              --thread {threads} \
              --qualified_quality_phred {params.q} \
              --unqualified_percent_limit {params.unqualified} \
              --length_required {params.min_length} \
              --correction {params.detect} > {log:q} 2>&1
            """
else:
    rule fastp:
        input:
            r1=lambda wc: SAMPLES[wc.sample]["read1"]
        output:
            r1=TRIM_R1_PATTERN,
            json=FASTP_JSON_PATTERN,
            html=FASTP_HTML_PATTERN
        log:
            f"{LOG_DIR}/fastp/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/fastp/{{sample}}.tsv"
        threads: 16
        resources:
            mem_mb=16000,
            runtime=240
        params:
            q=config.get("qc", {}).get("fastp", {}).get("qualified_quality_phred", 20),
            unqualified=config.get("qc", {}).get("fastp", {}).get("unqualified_percent_limit", 40),
            min_length=config.get("qc", {}).get("fastp", {}).get("length_required", 36)
        conda:
            QC_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"
            fastp \
              --in1 {input.r1:q} --out1 {output.r1:q} \
              --json {output.json:q} --html {output.html:q} \
              --thread {threads} \
              --qualified_quality_phred {params.q} \
              --unqualified_percent_limit {params.unqualified} \
              --length_required {params.min_length} > {log:q} 2>&1
            """


rule fastqc_trimmed:
    input:
        trimmed_read
    output:
        html=f"{RESULTS_DIR}/qc/fastqc/trimmed/{{sample}}.{{read}}_fastqc.html",
        zip=f"{RESULTS_DIR}/qc/fastqc/trimmed/{{sample}}.{{read}}_fastqc.zip"
    log:
        f"{LOG_DIR}/fastqc/trimmed/{{sample}}.{{read}}.log"
    benchmark:
        f"{BENCHMARK_DIR}/fastqc/trimmed/{{sample}}.{{read}}.tsv"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        QC_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.html:q})" "$(dirname {log:q})" "{resources.tmpdir}"
        tmp_dir=$(mktemp -d "{resources.tmpdir}/fastqc.trimmed.XXXXXX")
        trap 'rm -rf "$tmp_dir"' EXIT
        fastqc --threads {threads} --outdir "$tmp_dir" {input:q} > {log:q} 2>&1
        mv "$tmp_dir"/*_fastqc.html {output.html:q}
        mv "$tmp_dir"/*_fastqc.zip {output.zip:q}
        """


if STATE["screen_enabled"]:
    rule fastq_screen:
        input:
            r1=TRIM_R1_PATTERN
        output:
            txt=f"{RESULTS_DIR}/qc/fastq_screen/{{sample}}_screen.txt",
            html=f"{RESULTS_DIR}/qc/fastq_screen/{{sample}}_screen.html"
        log:
            f"{LOG_DIR}/fastq_screen/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/fastq_screen/{{sample}}.tsv"
        threads: 8
        resources:
            mem_mb=16000,
            runtime=240
        params:
            conf=STATE["screen_config"],
            subset=config.get("qc", {}).get("fastq_screen", {}).get("subset", 100000)
        conda:
            QC_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output.txt:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/fastq_screen.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            fastq_screen --conf {params.conf:q} --subset {params.subset} \
              --threads {threads} --outdir "$tmp_dir" {input.r1:q} > {log:q} 2>&1
            mv "$tmp_dir"/*_screen.txt {output.txt:q}
            mv "$tmp_dir"/*_screen.html {output.html:q}
            """


if MODE in {"bulk_rna", "dna_reseq"}:
    rule stage_reference_fasta:
        input:
            STATE["fasta"]
        output:
            REF_FASTA
        log:
            f"{LOG_DIR}/reference/stage_fasta.log"
        resources:
            mem_mb=1000,
            runtime=10
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
            ln -sfn "$(readlink -f {input:q})" {output:q}
            printf 'Staged %s -> %s\n' {input:q} {output:q} > {log:q}
            """

    rule reference_faidx:
        input:
            REF_FASTA
        output:
            f"{REF_FASTA}.fai"
        log:
            f"{LOG_DIR}/reference/faidx.log"
        threads: 2
        resources:
            mem_mb=4000,
            runtime=60
        conda:
            DNA_ENV
        shell:
            "samtools faidx {input:q} > {log:q} 2>&1"


if MODE == "bulk_rna":
    rule stage_reference_gtf:
        input:
            STATE["gtf"]
        output:
            REF_GTF
        log:
            f"{LOG_DIR}/reference/stage_gtf.log"
        resources:
            mem_mb=1000,
            runtime=10
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
            ln -sfn "$(readlink -f {input:q})" {output:q}
            printf 'Staged %s -> %s\n' {input:q} {output:q} > {log:q}
            """


rule multiqc:
    input:
        COMMON_QC_TARGETS + MODE_QC_TARGETS
    output:
        MULTIQC_REPORT
    log:
        f"{LOG_DIR}/multiqc/multiqc.log"
    benchmark:
        f"{BENCHMARK_DIR}/multiqc/multiqc.tsv"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=180
    conda:
        QC_ENV
    shell:
        r"""
        set -euo pipefail
        out_dir=$(dirname {output:q})
        mkdir -p "$out_dir" "$(dirname {log:q})" "{resources.tmpdir}"
        file_list=$(mktemp "{resources.tmpdir}/multiqc.inputs.XXXXXX")
        trap 'rm -f "$file_list"' EXIT
        printf '%s\n' {input:q} > "$file_list"
        multiqc --force --file-list "$file_list" --outdir "$out_dir" \
          --filename "$(basename {output:q})" > {log:q} 2>&1
        """


rule build_portal_manifest:
    input:
        multiqc=MULTIQC_REPORT,
        mode=MODE_TARGETS,
        fastp=FASTP_TARGETS
    output:
        manifest=PORTAL_MANIFEST,
        samples=f"{PORTAL_DIR}/tables/samples.tsv",
        qc=f"{PORTAL_DIR}/tables/sample_qc.tsv",
        mapping=f"{PORTAL_DIR}/tables/mapping_qc.tsv",
        config_used=f"{PROJECT_DIR}/provenance/config.used.yaml",
        checksums=f"{PROJECT_DIR}/provenance/result_checksums.sha256"
    log:
        f"{LOG_DIR}/portal/build_manifest.log"
    resources:
        mem_mb=8000,
        runtime=120
    params:
        project_dir=PROJECT_DIR,
        project_id=STATE["project_id"],
        display_name=STATE["display_name"],
        mode=MODE,
        pipeline_version=PIPELINE_VERSION,
        sample_sheet=STATE["sample_sheet"],
        contrast_ids=STATE["contrast_ids"],
        contrasts=STATE["contrasts"],
        config=config
    conda:
        REPORT_ENV
    script:
        f"{SCRIPT_DIR}/build_portal_manifest.py"
