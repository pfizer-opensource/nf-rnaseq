/*
Nextflow module to recombine MAJIQ outputs with metadata
*/

params.allow_missing_samples = true

process CREATE_CH_MAJIQ {

    executor 'local'
    
    input:
    val(nextflow_csv_row)
    val(majiq_out)
    // majiq_sj

    output:
    tuple val(majiq_group), val(meta), val(curr_majiq_sample_files), emit: majiq_out_meta
    
    exec:
    majiq_group = nextflow_csv_row.majiq_group
    meta = [:]
    meta.id = nextflow_csv_row.sample_name
    meta.majiq_group = nextflow_csv_row.majiq_group
    meta.strandedness = nextflow_csv_row.strandedness
    meta.orientation = nextflow_csv_row.orientation

    curr_majiq_sample_files = majiq_out.findAll {it =~ nextflow_csv_row.sample_name}
    
    if( curr_majiq_sample_files.isEmpty() ) {
        if(params.allow_missing_samples) {
            log.warn("No matching MAJIQ files found for ${nextflow_csv_row.sample_name}")
            curr_majiq_sample_files = null
            } else {
                exit 1, "No matching MAJIQ files found for ${nextflow_csv_row.sample_name}"
            }
    } else {
       curr_majiq_sample_files = curr_majiq_sample_files.sort()
    }
    
    meta.majiq_out = curr_majiq_sample_files
}