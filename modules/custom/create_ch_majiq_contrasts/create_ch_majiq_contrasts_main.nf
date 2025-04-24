/*
Nextflow module to recombine MAJIQ contrasts with metadata including MAJIQ groups and contrast identifiers

*/

params.allow_missing_samples = true

process CREATE_CH_MAJIQ_CONTRASTS {

    executor 'local'
    
    input:
    val(majiq_contrasts_row)
    val(nextflow_csv)
    val(meta)
    // majiq_sj

    output:
    tuple val(contrast_info), val(majiq_reference_samples), val(majiq_experiment_samples), emit: majiq_contrast_samples
    
    exec:
    contrast_info = [:]
    contrast_info.name = majiq_contrasts_row.contrast
    contrast_info.reference = majiq_contrasts_row.reference
    contrast_info.experiment = majiq_contrasts_row.experiment

    log.debug("Working on contrast ${majiq_contrasts_row.contrast} with reference condition: ${majiq_contrasts_row.reference}, experimental condition: ${majiq_contrasts_row.experiment}")
    
    // REFERENCE CONDITION //
    
    majiq_reference_meta = meta.findAll {it.majiq_group == majiq_contrasts_row.reference}
    
    if( majiq_reference_meta.isEmpty() ) {
        if(params.allow_missing_samples) {
            log.warn("No matching MAJIQ files found for ${majiq_contrasts_row.reference}")
            majiq_reference_samples = null
            } else {
                exit 1, "No matching MAJIQ files found for ${majiq_contrasts_row.reference}"
            }
    } else {
         majiq_reference_samples = majiq_reference_meta.majiq_out
    }
    log.debug("Found reference ${majiq_reference_samples}")

    
    
    // EXPERIMENTAL CONDITION //
    
    majiq_experiment_meta = meta.findAll {it.majiq_group == majiq_contrasts_row.experiment}
    
    if( majiq_experiment_meta.isEmpty() ) {
        if(params.allow_missing_samples) {
            log.warn("No matching MAJIQ files found for ${majiq_contrasts_row.experiment}")
            majiq_experiment_samples = null
            } else {
                exit 1, "No matching MAJIQ files found for ${majiq_contrasts_row.experiment}"
            }
    } else {
          majiq_experiment_samples = majiq_experiment_meta.majiq_out
    }
    log.debug("Found experiment ${majiq_experiment_samples}")
}
