/*
Nextlow module to run MultiQC

Requires: none

Params: None

Input:
  -- one path entry for every tool collected in the pipeline using a named input such as `path('<tool name>/*')`
     to prevent collisions; this must be modified to match the needs of each pipeline this module is used in

Output:
  -- path to the multiqc html report, emitted as `report`
  -- path to all multiqc *_data directories, emitted as `data`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

nextflow.enable.dsl = 2

process MULTIQC {
    label 'process_1cpu_1gb_30min'

    input:
    path multiqc_config
    path software_versions
    path('fastqc_raw/*')
    path('fastp/*')
    path('fastqc_postfastp/*')
    path('star/*')
    path('salmon/*')
    path('picard/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    """
    multiqc -f $args .

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        multiqc: \$(echo \$(multiqc --version) | sed -e 's/multiqc//g' -e 's/version//g' -e 's/,//g' -e 's/ \\+//g')
    END_VERSION
    """
}
