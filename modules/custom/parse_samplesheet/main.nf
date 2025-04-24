/*
Custom Nextflow module to parse an Illumina format SampleSheet.csv

Input:
  -- csv: an Illumina SampleSheet.csv

Output:
  -- SampleSheetNextflow.csv: the parsed sample sheet formatted with one sample per row
       to facilitate use with other Nextflow pipeline modules, emitted as `nextflow_csv`
  -- _version.yaml containing software version information gathered at run-time, emitted as `version`
*/

process PARSE_SAMPLESHEET {
    input:
    file csv

    output:
    path 'SampleSheetNextflow.csv', emit: nextflow_csv
    path '_version.yaml', emit: version

    """
    ${workflow.projectDir}/bin/parse_sample_sheet.py --samplesheet $csv

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        parse_sample_sheet.py: \$(echo \$(parse_sample_sheet.py --version) | sed 's/parse_sample_sheet.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}   
