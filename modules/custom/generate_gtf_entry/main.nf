/*
Nextflow module to create GTF rows from a YAML specification

Requires: generate_gtf_entry.py

Params: none

Input:
  -- yaml, path to a YAML file specifying custom gene structure

Output:
  -- tuple of metadata and path to the generated GTF file, emitted as `gtf`
  -- _versions.yaml file capturing the tool version, emitted as `version`

Customization: publishDir in modules.config
*/

nextflow.enable.dsl=2

process GENERATE_GTF_ENTRY {
    executor 'local'

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path('*.gtf'), emit: gtf
    path '_version.yaml', emit: version

    script:
    """
    /home/moccir/projects/rdru-nextflow/bin/generate_gtf_entry.py $yaml > custom.gtf

    cat <<-END_VERSION >_version.yaml
    "${task.process}":
        generate_gtf_entry.py: \$(/home/moccir/projects/rdru-nextflow/bin/generate_gtf_entry.py --version 2>&1 | sed 's/generate_gtf_entry.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
