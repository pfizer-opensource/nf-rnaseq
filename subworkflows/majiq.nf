nextflow.enable.dsl = 2

// parameters
// allow_missing_samples: do not exit with error if majiq files are not
//      found for all samples in SampleSheet

params.allow_missing_samples = true
params.majiq_build_dir = "./processed_nextflow_majiq/majiq/build_out"

include { GET_MAJIQ_SETTINGS } from '../modules/custom/get_majiq_settings/get_majiq_settings_main'
include { MAJIQ_BUILD } from '../modules/majiq/build/main'
include { CREATE_CH_MAJIQ } from '../modules/custom/create_ch_majiq/create_ch_majiq_main'
include { MAJIQ_SAMPLE_PSI } from '../modules/majiq/sample_psi/sample_psi_main'
include { CREATE_CH_MAJIQ_CONTRASTS } from '../modules/custom/create_ch_majiq_contrasts/create_ch_majiq_contrasts_main'
include { MAJIQ_GROUP_PSI } from '../modules/majiq/group_psi/group_psi_main'
include { MAJIQ_DELTAPSI } from '../modules/majiq/deltapsi/deltapsi_main'
// include { VOILA_TSV } from '../modules/majiq/voila_tsv/voila_tsv_main'

workflow MAJIQ {

    take:
    nextflow_csv
    majiq_contrasts
    genome_identifier
    bam_sorted
    bam_index
    gff3
    
    main:
    
    ch_versions = Channel.empty()
    
    GET_MAJIQ_SETTINGS( nextflow_csv, majiq_contrasts, genome_identifier.tokenize( ':' ).map{ it[0] })
    ch_versions = ch_versions.mix(GET_MAJIQ_SETTINGS.out.version.first())
    
    MAJIQ_BUILD( GET_MAJIQ_SETTINGS.out.majiq_settings, bam_sorted, bam_index, gff3)
    ch_versions = ch_versions.mix(MAJIQ_BUILD.out.version.first())
    
    CREATE_CH_MAJIQ( 
        nextflow_csv.splitCsv(header: true).map{ it }, 
        MAJIQ_BUILD.out.majiq_out.first()
        )  
    CREATE_CH_MAJIQ.out.majiq_out_meta.map({it -> [it[0], it[1], it[2]]}).set{ CH_MAJIQ_GROUPNAMES }
    CH_MAJIQ_GROUPNAMES.groupTuple(by: 0).set{ CH_MAJIQ_GROUPS }

    MAJIQ_SAMPLE_PSI(MAJIQ_BUILD.out.majiq_out.flatten())
    ch_versions = ch_versions.mix(MAJIQ_SAMPLE_PSI.out.version.first())
    
    // MAJIQ_GROUP_PSI(CH_MAJIQ_GROUPS.map{ it[0] },CH_MAJIQ_GROUPS.map{ it[2] }.flatten().collect())
    // ch_versions = ch_versions.mix(MAJIQ_GROUP_PSI.out.version.first())
     
    CREATE_CH_MAJIQ_CONTRASTS( 
       majiq_contrasts.splitCsv(header: true).map{ it }, 
       nextflow_csv.splitCsv(header: true).collect(), 
       CREATE_CH_MAJIQ.out.majiq_out_meta.map{ it[1]}.collect()
       )
    CREATE_CH_MAJIQ_CONTRASTS.out.majiq_contrast_samples.set{ CH_MAJIQ_CONTRASTS }
    
    // CH_MAJIQ_CONTRASTS.view()
   
    MAJIQ_DELTAPSI(CH_MAJIQ_CONTRASTS.map{ it[0] },CH_MAJIQ_CONTRASTS.map{ it[1].flatten() }, CH_MAJIQ_CONTRASTS.map{ it[2].flatten() },MAJIQ_BUILD.out.splicegraph.first())
    ch_versions = ch_versions.mix(MAJIQ_DELTAPSI.out.version.first())
    
    // VOILA_TSV(MAJIQ_BUILD.out.splicegraph.first(), MAJIQ_DELTAPSI.out.deltapsi_voila)
    // ch_versions = ch_versions.mix(VOILA_TSV.out.version.first())
    
    emit:
    versions = ch_versions
}