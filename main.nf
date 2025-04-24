#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Validate pipeline parameters
WorkflowMain.initialize(workflow, params, log)

// Perform additional pipeline-specific inputs
WorkflowRnaseq.initialize(params, log)

include { PARSE_SAMPLESHEET } from './modules/custom/parse_samplesheet/main'
include { create_fastq_entry } from './bin/utilities.groovy'
include { LOAD_GENOME } from './subworkflows/load_genome'
include { FASTQC } from './modules/fastqc/main'
include { FASTP } from './modules/fastp/main'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc/main'
include { STAR_ALIGN } from './modules/star/align/main'
include { STAR_2PASS } from './subworkflows/star_2pass'
include { SAMTOOLS_INDEX } from './modules/samtools/index/main'
include { PICARD_COLLECTRNASEQMETRICS } from './modules/picard/collectrnaseqmetrics/main'
include { SALMON_QUANT_ALIGNED } from './modules/salmon/quant_aligned/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { MAKE_EXPERIMENT_JSON } from './modules/custom/make_experiment_json/main'
include { MULTIQC } from './modules/multiqc/main'
include { MAJIQ } from './subworkflows/majiq'


// create channels
ch_samplesheet = Channel
    .fromPath( params.input + '/SampleSheet.csv', checkIfExists: true )


Channel
    .value(params.junctions)
    .set{ ch_junctions }
    
Channel
    .value(params.genome_id)
    .set{ ch_genome_identifier }
    
if(params.run_majiq) {
    Channel
        .fromPath(params.input + '/majiq_contrasts.csv', checkIfExists: true )
        .set{ ch_majiq_contrasts }
    }
    
if( params.junctions != 'NO_FILE' ) {
    ch_junctions = Channel
        .fromPath(params.junctions, checkIfExists: true)
        .collect()
} else {
    ch_junctions = file(params.junctions)
}

ch_versions = Channel.empty()

ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yaml", checkIfExists: true)

workflow {
    PARSE_SAMPLESHEET( ch_samplesheet )
    ch_versions = ch_versions.mix(PARSE_SAMPLESHEET.out.version)
    ch_reads = PARSE_SAMPLESHEET.out.nextflow_csv
        .splitCsv(header: true)
        .map{ create_fastq_entry(it) }

    LOAD_GENOME()
    ch_versions = ch_versions.mix(LOAD_GENOME.out.versions)

    FASTQC( ch_reads )
    ch_versions = ch_versions.mix( FASTQC.out.version )

    if (params.trim_reads) {
      FASTP( ch_reads )
      ch_versions = ch_versions.mix( FASTP.out.version )
      ch_fastpout = FASTP.out.fastp_json.collect{ it[1] }
      FASTQC_TRIMMED( FASTP.out.fastp_reads )
      ch_alignable_reads = FASTP.out.fastp_reads
      ch_versions = ch_versions.mix( FASTQC_TRIMMED.out.version )
      ch_fastqc_trimmed = FASTQC_TRIMMED.out.fastqc_dirs.collect{ it[1] }
    } else {
        ch_alignable_reads = ch_reads
        Channel.value( [] ).set{ ch_fastpout }
        Channel.value( [] ).set{ ch_fastqc_trimmed }
    }

    if (params.star_2pass) {
        STAR_2PASS(
            ch_alignable_reads,
            LOAD_GENOME.out.star_index.map{ it -> it[1] },
            LOAD_GENOME.out.gtf.map{ it -> it[1] },
            ch_junctions,
            LOAD_GENOME.out.genome_fasta.map{ it -> it[1] } )
        ch_bam_transcriptome = STAR_2PASS.out.bam_transcriptome
        ch_star_logs = STAR_2PASS.out.log_final.collect{ it[1] }
        ch_versions = ch_versions.mix( STAR_2PASS.out.versions)

        SAMTOOLS_INDEX( STAR_2PASS.out.bam_sorted )
        ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.version )

        ch_bam = STAR_2PASS.out.bam_sorted_pass1
            .map{ it ->
                    new_meta = it[0].clone()
                    new_meta.id = "${new_meta.id}.pass1"
                    tuple(new_meta, it[1])
                    }
            .mix( STAR_2PASS.out.bam_sorted )

    } else {
        STAR_ALIGN(
            ch_alignable_reads,
            LOAD_GENOME.out.star_index.map{ it -> it[1] },
            LOAD_GENOME.out.gtf.map{ it -> it[1] },
            ch_junctions
        )
        ch_bam = STAR_ALIGN.out.bam_sorted
        ch_bam_transcriptome = STAR_ALIGN.out.bam_transcriptome
        ch_star_logs = STAR_ALIGN.out.log_final.collect{ it[1] }
        ch_versions = ch_versions.mix(STAR_ALIGN.out.version)

        SAMTOOLS_INDEX( STAR_ALIGN.out.bam_sorted )
        ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.version )
    }

    PICARD_COLLECTRNASEQMETRICS(
        ch_bam,
        LOAD_GENOME.out.refflat.map{ it -> it[1] },
        LOAD_GENOME.out.rrna_intervals.map{ it -> it[1] }
        )
    ch_versions = ch_versions.mix( PICARD_COLLECTRNASEQMETRICS.out.version )

    SALMON_QUANT_ALIGNED( ch_bam_transcriptome, LOAD_GENOME.out.transcriptome_fasta.map{ it -> it[1] } )
    ch_versions = ch_versions.mix( SALMON_QUANT_ALIGNED.out.version )

//    SAMTOOLS_INDEX.out.indexed_alignment.collect{ it[0] }.view()

    if(params.run_majiq) {
        MAJIQ(
            PARSE_SAMPLESHEET.out.nextflow_csv,
            ch_majiq_contrasts,
            ch_genome_identifier,
	          SAMTOOLS_INDEX.out.indexed_alignment.collect{ it[1] },
	          SAMTOOLS_INDEX.out.indexed_alignment.collect{ it[2] },
	          LOAD_GENOME.out.gff3.map{ it -> it[1] }
	          )
        ch_versions = ch_versions.mix(MAJIQ.out.versions)
        }


    CUSTOM_DUMPSOFTWAREVERSIONS( ch_versions.unique().collectFile(name: 'collated_versions.yml') )

    MAKE_EXPERIMENT_JSON( ch_reads.collect(flat: false), LOAD_GENOME.out.genome_fasta )
    
    
    MULTIQC( 
        ch_multiqc_config,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        FASTQC.out.fastqc_dirs.collect{ it[1] },
        ch_fastpout,
        ch_fastqc_trimmed,
        ch_star_logs,
        SALMON_QUANT_ALIGNED.out.salmon_dirs.collect{ it[1] },
        PICARD_COLLECTRNASEQMETRICS.out.picard_rna_metrics.collect{ it[1] }
        )
}
