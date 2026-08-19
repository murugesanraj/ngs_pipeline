"""Short-read de novo genome and transcriptome assembly branches."""


ASSEMBLER = str(config.get("assembly", {}).get("genome_assembler", "spades"))
MIN_CONTIG = int(config.get("assembly", {}).get("minimum_contig_length", 500))
if MODE == "denovo_genome" and ASSEMBLER not in {"spades", "megahit"}:
    raise WorkflowError("assembly.genome_assembler must be spades or megahit")


if MODE == "denovo_genome" and ASSEMBLER == "spades" and PAIRED:
    rule assemble_spades:
        input:
            r1=TRIM_R1_PATTERN,
            r2=TRIM_R2_PATTERN
        output:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        log:
            f"{LOG_DIR}/assembly/{{sample}}.spades.log"
        benchmark:
            f"{BENCHMARK_DIR}/assembly/{{sample}}.spades.tsv"
        threads: 32
        resources:
            mem_mb=240000,
            runtime=2880
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/spades.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            memory_gb=$(( {resources.mem_mb} / 1024 ))
            spades.py --careful --threads {threads} --memory "$memory_gb" \
              -1 {input.r1:q} -2 {input.r2:q} -o "$tmp_dir/out" > {log:q} 2>&1
            seqkit seq -m {MIN_CONTIG} "$tmp_dir/out/scaffolds.fasta" > {output:q}
            test -s {output:q}
            """

elif MODE == "denovo_genome" and ASSEMBLER == "spades":
    rule assemble_spades:
        input:
            r1=TRIM_R1_PATTERN
        output:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        log:
            f"{LOG_DIR}/assembly/{{sample}}.spades.log"
        benchmark:
            f"{BENCHMARK_DIR}/assembly/{{sample}}.spades.tsv"
        threads: 32
        resources:
            mem_mb=240000,
            runtime=2880
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/spades.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            memory_gb=$(( {resources.mem_mb} / 1024 ))
            spades.py --careful --threads {threads} --memory "$memory_gb" \
              -s {input.r1:q} -o "$tmp_dir/out" > {log:q} 2>&1
            seqkit seq -m {MIN_CONTIG} "$tmp_dir/out/scaffolds.fasta" > {output:q}
            test -s {output:q}
            """

elif MODE == "denovo_genome" and ASSEMBLER == "megahit" and PAIRED:
    rule assemble_megahit:
        input:
            r1=TRIM_R1_PATTERN,
            r2=TRIM_R2_PATTERN
        output:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        log:
            f"{LOG_DIR}/assembly/{{sample}}.megahit.log"
        benchmark:
            f"{BENCHMARK_DIR}/assembly/{{sample}}.megahit.tsv"
        threads: 32
        resources:
            mem_mb=240000,
            runtime=2880
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/megahit.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            megahit -1 {input.r1:q} -2 {input.r2:q} -t {threads} \
              --min-contig-len {MIN_CONTIG} -o "$tmp_dir/out" > {log:q} 2>&1
            cp "$tmp_dir/out/final.contigs.fa" {output:q}
            test -s {output:q}
            """

elif MODE == "denovo_genome" and ASSEMBLER == "megahit":
    rule assemble_megahit:
        input:
            r1=TRIM_R1_PATTERN
        output:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        log:
            f"{LOG_DIR}/assembly/{{sample}}.megahit.log"
        benchmark:
            f"{BENCHMARK_DIR}/assembly/{{sample}}.megahit.tsv"
        threads: 32
        resources:
            mem_mb=240000,
            runtime=2880
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/megahit.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            megahit -r {input.r1:q} -t {threads} --min-contig-len {MIN_CONTIG} \
              -o "$tmp_dir/out" > {log:q} 2>&1
            cp "$tmp_dir/out/final.contigs.fa" {output:q}
            test -s {output:q}
            """

elif MODE == "denovo_transcriptome" and PAIRED:
    rule assemble_trinity:
        input:
            r1=TRIM_R1_PATTERN,
            r2=TRIM_R2_PATTERN
        output:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        log:
            f"{LOG_DIR}/assembly/{{sample}}.trinity.log"
        benchmark:
            f"{BENCHMARK_DIR}/assembly/{{sample}}.trinity.tsv"
        threads: 32
        resources:
            mem_mb=240000,
            runtime=2880
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/trinity.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            memory_gb=$(( {resources.mem_mb} / 1024 ))
            Trinity --seqType fq --CPU {threads} --max_memory "${{memory_gb}}G" \
              --left {input.r1:q} --right {input.r2:q} \
              --output "$tmp_dir/trinity_out" > {log:q} 2>&1
            seqkit seq -m {MIN_CONTIG} "$tmp_dir/trinity_out.Trinity.fasta" > {output:q}
            test -s {output:q}
            """

elif MODE == "denovo_transcriptome":
    rule assemble_trinity:
        input:
            r1=TRIM_R1_PATTERN
        output:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        log:
            f"{LOG_DIR}/assembly/{{sample}}.trinity.log"
        benchmark:
            f"{BENCHMARK_DIR}/assembly/{{sample}}.trinity.tsv"
        threads: 32
        resources:
            mem_mb=240000,
            runtime=2880
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/trinity.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            memory_gb=$(( {resources.mem_mb} / 1024 ))
            Trinity --seqType fq --CPU {threads} --max_memory "${{memory_gb}}G" \
              --single {input.r1:q} --output "$tmp_dir/trinity_out" > {log:q} 2>&1
            seqkit seq -m {MIN_CONTIG} "$tmp_dir/trinity_out.Trinity.fasta" > {output:q}
            test -s {output:q}
            """


rule quast:
    input:
        f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
    output:
        report=f"{RESULTS_DIR}/assembly/{{sample}}/quast/report.tsv",
        html=f"{RESULTS_DIR}/assembly/{{sample}}/quast/report.html"
    log:
        f"{LOG_DIR}/quast/{{sample}}.log"
    benchmark:
        f"{BENCHMARK_DIR}/quast/{{sample}}.tsv"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=480
    params:
        out_dir=lambda wc: f"{RESULTS_DIR}/assembly/{wc.sample}/quast",
        mode="--rna-finding" if MODE == "denovo_transcriptome" else ""
    conda:
        ASSEMBLY_ENV
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.out_dir:q} "$(dirname {log:q})"
        quast.py --threads {threads} {params.mode} --output-dir {params.out_dir:q} \
          {input:q} > {log:q} 2>&1
        test -s {output.report:q}
        test -s {output.html:q}
        """


if STATE["busco_enabled"]:
    rule busco:
        input:
            f"{RESULTS_DIR}/assembly/{{sample}}/contigs.fasta"
        output:
            touch(f"{RESULTS_DIR}/assembly/{{sample}}/busco/.complete")
        log:
            f"{LOG_DIR}/busco/{{sample}}.log"
        benchmark:
            f"{BENCHMARK_DIR}/busco/{{sample}}.tsv"
        threads: 16
        resources:
            mem_mb=64000,
            runtime=1440
        params:
            lineage=STATE["busco_lineage"],
            mode="transcriptome" if MODE == "denovo_transcriptome" else "genome",
            out_dir=lambda wc: f"{RESULTS_DIR}/assembly/{wc.sample}/busco"
        conda:
            ASSEMBLY_ENV
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {params.out_dir:q})" "$(dirname {log:q})" "{resources.tmpdir}"
            tmp_dir=$(mktemp -d "{resources.tmpdir}/busco.{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$tmp_dir"' EXIT
            cd "$tmp_dir"
            busco --offline --cpu {threads} --mode {params.mode:q} \
              --lineage_dataset {params.lineage:q} --in {input:q} --out busco > {log:q} 2>&1
            rm -rf {params.out_dir:q}
            mv "$tmp_dir/busco" {params.out_dir:q}
            touch {output:q}
            """
