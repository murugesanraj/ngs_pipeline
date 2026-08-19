/*
 * Hellbender NGS Pipeline -- Nextflow DSL2 entry point.
 *
 * Project settings are supplied with `-params-file project/config.yaml`.
 */

params {
    pipeline: Map
    project: Map
    analysis: Map
    reference: Map = [:]
    qc: Map = [:]
    rna: Map = [:]
    dna: Map = [:]
    assembly: Map = [:]
    portal: Map = [:]
}

include {
    FASTQC_RAW
    FASTP_PE
    FASTP_SE
    FASTQC_TRIMMED
    FASTQ_SCREEN
    MULTIQC
    BUILD_PORTAL
} from './nextflow/modules/common'

include {
    STAR_INDEX
    STAR_ALIGN_PE
    STAR_ALIGN_SE
    RNA_BIGWIG
    FEATURECOUNTS
    CLEAN_FEATURECOUNTS
    DESEQ2
} from './nextflow/modules/rna'

include {
    REFERENCE_FAIDX
    BWA_INDEX
    BWA_INDEX_PREBUILT
    BWA_ALIGN_MARKDUP_PE
    BWA_ALIGN_MARKDUP_SE
    MOSDEPTH
    CALL_VARIANTS
} from './nextflow/modules/dna'

include {
    ASSEMBLE_SPADES_PE
    ASSEMBLE_SPADES_SE
    ASSEMBLE_MEGAHIT_PE
    ASSEMBLE_MEGAHIT_SE
    ASSEMBLE_TRINITY_PE
    ASSEMBLE_TRINITY_SE
    QUAST
    BUSCO
} from './nextflow/modules/assembly'


def requireChoice(value, choices, label) {
    def normalized = value == null ? '' : value.toString()
    if (!choices.contains(normalized)) {
        error "${label} must be one of: ${choices.join(', ')}"
    }
    normalized
}


def mapValue(values, key, fallback) {
    values != null && values.containsKey(key) && values[key] != null ? values[key] : fallback
}


def boolValue(value) {
    if (value instanceof Boolean) {
        return value
    }
    value.toString().toBoolean()
}


workflow {
    def modes = ['bulk_rna', 'dna_reseq', 'denovo_genome', 'denovo_transcriptome']
    def layouts = ['paired', 'single']
    def mode = requireChoice(params.analysis.mode, modes, 'analysis.mode')
    def layout = requireChoice(params.analysis.layout, layouts, 'analysis.layout')
    def paired = layout == 'paired'
    def outdir = file(params.project.output_dir).toAbsolutePath().toString()
    def projectId = params.project.id.toString()

    if (!(projectId ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,99}/)) {
        error 'project.id contains unsafe characters'
    }

    sample_sheet_file = file(params.project.sample_sheet, checkIfExists: true)
    sample_rows = channel.fromPath(params.project.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)

    if (paired) {
        samples = sample_rows.map { row ->
            def sampleId = row.sample_id == null ? '' : row.sample_id.toString()
            if (!(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,99}/)) {
                error "Unsafe or missing sample_id: ${sampleId}"
            }
            if (row.read1 == null || row.read1.toString().isEmpty()) {
                error "Sample ${sampleId} has no read1"
            }
            if (row.read2 == null || row.read2.toString().isEmpty()) {
                error "Paired sample ${sampleId} has no read2"
            }
            def meta = row.findAll { key, _value -> key != 'read1' && key != 'read2' }
            tuple(
                meta,
                file(row.read1.toString(), checkIfExists: true),
                file(row.read2.toString(), checkIfExists: true)
            )
        }
        raw_reads = samples.flatMap { meta, read1, read2 ->
            [tuple(meta, 'R1', read1), tuple(meta, 'R2', read2)]
        }
    } else {
        samples = sample_rows.map { row ->
            def sampleId = row.sample_id == null ? '' : row.sample_id.toString()
            if (!(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,99}/)) {
                error "Unsafe or missing sample_id: ${sampleId}"
            }
            if (row.read1 == null || row.read1.toString().isEmpty()) {
                error "Sample ${sampleId} has no read1"
            }
            if (row.read2 != null && !row.read2.toString().isEmpty()) {
                error "Single-end sample ${sampleId} unexpectedly has read2"
            }
            def meta = row.findAll { key, _value -> key != 'read1' && key != 'read2' }
            tuple(meta, file(row.read1.toString(), checkIfExists: true))
        }
        raw_reads = samples.map { meta, read1 -> tuple(meta, 'R1', read1) }
    }

    def fastpConfig = mapValue(params.qc, 'fastp', [:])
    def quality = mapValue(fastpConfig, 'qualified_quality_phred', 20) as Integer
    def unqualified = mapValue(fastpConfig, 'unqualified_percent_limit', 40) as Integer
    def minLength = mapValue(fastpConfig, 'length_required', 36) as Integer
    def detectAdapters = boolValue(mapValue(fastpConfig, 'detect_adapter_for_pe', true))

    raw_qc = FASTQC_RAW(raw_reads, outdir)
    if (paired) {
        fastp = FASTP_PE(samples, outdir, quality, unqualified, minLength, detectAdapters)
        trimmed_reads = fastp.reads.flatMap { meta, read1, read2 ->
            [tuple(meta, 'R1', read1), tuple(meta, 'R2', read2)]
        }
        screen_reads = fastp.reads.map { meta, read1, _read2 -> tuple(meta, read1) }
    } else {
        fastp = FASTP_SE(samples, outdir, quality, unqualified, minLength)
        trimmed_reads = fastp.reads.map { meta, read1 -> tuple(meta, 'R1', read1) }
        screen_reads = fastp.reads.map { meta, read1 -> tuple(meta, read1) }
    }
    trimmed_qc = FASTQC_TRIMMED(trimmed_reads, outdir)

    def screenConfig = mapValue(params.qc, 'fastq_screen', [:])
    def screenEnabled = boolValue(mapValue(screenConfig, 'enabled', false))
    if (screenEnabled) {
        screen_config_file = file(screenConfig.config, checkIfExists: true)
        screen = FASTQ_SCREEN(
            screen_reads,
            screen_config_file,
            outdir,
            mapValue(screenConfig, 'subset', 100000) as Integer
        )
        screen_qc_files = screen.reports.map { _meta, textReport, _htmlReport -> textReport }
    } else {
        screen_qc_files = channel.empty()
    }

    raw_qc_files = raw_qc.reports.map { _meta, _readLabel, _htmlReport, zipReport -> zipReport }
    fastp_qc_files = fastp.reports.map { _meta, jsonReport, _htmlReport -> jsonReport }
    trimmed_qc_files = trimmed_qc.reports.map { _meta, _readLabel, _htmlReport, zipReport -> zipReport }
    common_qc_files = raw_qc_files.mix(fastp_qc_files, trimmed_qc_files, screen_qc_files)

    clean_counts_script = file(
        "${projectDir}/workflow/scripts/aggregate_featurecounts_nextflow.py",
        checkIfExists: true
    )
    deseq2_script = file(
        "${projectDir}/workflow/scripts/deseq2_nextflow.R",
        checkIfExists: true
    )
    portal_script = file(
        "${projectDir}/workflow/scripts/build_portal_manifest_nextflow.py",
        checkIfExists: true
    )

    if (mode == 'bulk_rna') {
        fasta = file(params.reference.fasta, checkIfExists: true)
        gtf = file(params.reference.annotation_gtf, checkIfExists: true)
        def overhang = mapValue(params.rna, 'star_sjdb_overhang', 149) as Integer
        if (mapValue(params.reference, 'star_index', null)) {
            star_index = file(params.reference.star_index, checkIfExists: true)
        } else {
            built_star_index = STAR_INDEX(fasta, gtf, outdir, overhang)
            star_index = built_star_index.index
        }

        if (paired) {
            star = STAR_ALIGN_PE(fastp.reads, star_index, outdir)
        } else {
            star = STAR_ALIGN_SE(fastp.reads, star_index, outdir)
        }
        bigwig = RNA_BIGWIG(
            star.alignment.map { meta, bam, bai, _finalLog, _junctions, _geneCounts ->
                tuple(meta, bam, bai)
            },
            outdir
        )
        bam_files = star.alignment.map { _meta, bam, _bai, _finalLog, _junctions, _geneCounts -> bam }
            .collect()
        featurecounts = FEATURECOUNTS(
            bam_files,
            gtf,
            outdir,
            paired,
            mapValue(params.rna, 'strandedness', 'unstranded').toString(),
            mapValue(params.rna, 'featurecounts_feature_type', 'exon').toString(),
            mapValue(params.rna, 'featurecounts_id_attribute', 'gene_id').toString()
        )
        clean_counts = CLEAN_FEATURECOUNTS(
            featurecounts.raw,
            clean_counts_script,
            outdir
        )

        rna_logs = star.alignment.map { _meta, _bam, _bai, finalLog, _junctions, _geneCounts -> finalLog }
        rna_tracks = bigwig.tracks.map { _meta, track -> track }
        rna_analysis_base = rna_logs.mix(clean_counts.counts, rna_tracks)
        mode_qc_files = rna_logs.mix(featurecounts.summary)

        def contrastsSetting = mapValue(params.rna, 'contrasts', null)
        if (contrastsSetting != null && !contrastsSetting.toString().isEmpty()) {
            contrasts_file = file(contrastsSetting, checkIfExists: true)
            deseq2 = DESEQ2(
                clean_counts.counts,
                sample_sheet_file,
                contrasts_file,
                deseq2_script,
                outdir,
                mapValue(params.rna, 'design', ['condition'])
                    .collect { term -> term.toString() }
                    .join(','),
                mapValue(params.rna, 'alpha', 0.05) as Double,
                mapValue(params.rna, 'min_total_count', 10) as Integer
            )
            analysis_files = rna_analysis_base.mix(deseq2.results, channel.value(contrasts_file))
        } else {
            analysis_files = rna_analysis_base
        }
    } else if (mode == 'dna_reseq') {
        fasta = file(params.reference.fasta, checkIfExists: true)
        faidx = REFERENCE_FAIDX(fasta, outdir)
        def bwaPrefixSetting = mapValue(params.reference, 'bwa_index_prefix', null)
        if (bwaPrefixSetting != null && !bwaPrefixSetting.toString().isEmpty()) {
            def suffixes = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
            def indexPaths = suffixes.collect { suffix ->
                file("${bwaPrefixSetting}${suffix}", checkIfExists: true)
            }
            prebuilt_bwa = BWA_INDEX_PREBUILT(
                indexPaths,
                file(bwaPrefixSetting).name,
                fasta
            )
            bwa_index = prebuilt_bwa.index
        } else {
            built_bwa = BWA_INDEX(fasta, outdir)
            bwa_index = built_bwa.index
        }

        if (paired) {
            bwa = BWA_ALIGN_MARKDUP_PE(fastp.reads, bwa_index, outdir)
        } else {
            bwa = BWA_ALIGN_MARKDUP_SE(fastp.reads, bwa_index, outdir)
        }
        mosdepth = MOSDEPTH(
            bwa.alignment.map { meta, bam, bai, _flagstat, _idxstats, _duplicateStats ->
                tuple(meta, bam, bai)
            },
            outdir
        )
        flagstats = bwa.alignment.map { _meta, _bam, _bai, flagstat, _idxstats, _duplicateStats -> flagstat }
        coverage_summaries = mosdepth.coverage.map { _meta, summary, _distribution -> summary }
        dna_analysis_base = flagstats.mix(coverage_summaries)
        mode_qc_files = flagstats.mix(coverage_summaries)

        if (boolValue(mapValue(params.dna, 'call_variants', true))) {
            variants = CALL_VARIANTS(
                bwa.alignment.map { meta, bam, bai, _flagstat, _idxstats, _duplicateStats ->
                    tuple(meta, bam, bai)
                },
                fasta,
                faidx.index,
                outdir,
                mapValue(params.dna, 'minimum_variant_quality', 20) as Double,
                mapValue(params.dna, 'minimum_depth', 5) as Integer
            )
            variant_files = variants.results.flatMap { _meta, vcf, _index, stats -> [vcf, stats] }
            variant_stats = variants.results.map { _meta, _vcf, _index, stats -> stats }
            analysis_files = dna_analysis_base.mix(variant_files)
            mode_qc_files = mode_qc_files.mix(variant_stats)
        } else {
            analysis_files = dna_analysis_base
        }
    } else {
        def minContig = mapValue(params.assembly, 'minimum_contig_length', 500) as Integer
        if (mode == 'denovo_genome') {
            def assembler = requireChoice(
                mapValue(params.assembly, 'genome_assembler', 'spades'),
                ['spades', 'megahit'],
                'assembly.genome_assembler'
            )
            if (assembler == 'spades' && paired) {
                assembly_result = ASSEMBLE_SPADES_PE(fastp.reads, outdir, minContig)
            } else if (assembler == 'spades') {
                assembly_result = ASSEMBLE_SPADES_SE(fastp.reads, outdir, minContig)
            } else if (paired) {
                assembly_result = ASSEMBLE_MEGAHIT_PE(fastp.reads, outdir, minContig)
            } else {
                assembly_result = ASSEMBLE_MEGAHIT_SE(fastp.reads, outdir, minContig)
            }
        } else if (paired) {
            assembly_result = ASSEMBLE_TRINITY_PE(fastp.reads, outdir, minContig)
        } else {
            assembly_result = ASSEMBLE_TRINITY_SE(fastp.reads, outdir, minContig)
        }

        quast = QUAST(assembly_result.contigs, outdir, mode == 'denovo_transcriptome')
        contig_files = assembly_result.contigs.map { _meta, contigs -> contigs }
        quast_files = quast.reports.flatMap { _meta, tableReport, htmlReport ->
            [tableReport, htmlReport]
        }
        analysis_files = contig_files.mix(quast_files)
        mode_qc_files = quast.reports.map { _meta, tableReport, _htmlReport -> tableReport }

        def buscoConfig = mapValue(params.assembly, 'busco', [:])
        if (boolValue(mapValue(buscoConfig, 'enabled', false))) {
            lineage = file(buscoConfig.lineage_path, checkIfExists: true)
            busco = BUSCO(
                assembly_result.contigs,
                lineage,
                outdir,
                mode == 'denovo_transcriptome' ? 'transcriptome' : 'genome'
            )
            analysis_files = analysis_files.mix(busco.results.map { _meta, directory -> directory })
        }
    }

    all_qc_files = common_qc_files.mix(mode_qc_files)
    multiqc = MULTIQC(all_qc_files.collect(), outdir)
    fastp_jsons = fastp.reports.map { _meta, jsonReport, _htmlReport -> jsonReport }.collect()
    analysis_inputs = analysis_files.collect()
    config_json_text = groovy.json.JsonOutput.toJson([
        pipeline: params.pipeline,
        project: params.project,
        analysis: params.analysis,
        reference: params.reference,
        qc: params.qc,
        rna: params.rna,
        dna: params.dna,
        assembly: params.assembly,
        portal: params.portal
    ])
    config_json_file = channel.of(config_json_text)
        .collectFile(name: "${projectId}.config.json", newLine: false)
    BUILD_PORTAL(
        sample_sheet_file,
        config_json_file,
        mode,
        fastp_jsons,
        analysis_inputs,
        multiqc.report,
        portal_script,
        outdir
    )
}
