/*
Nextflow module to run Portcullis junctools set

Requires: none

Params: none

Input:
  -- tuple of metadata (Groovy map) and input_files (path), path to the files to be merged
  -- mode (string): set operation to apply [intersection, union, consensus, subtract, symmetric_difference,
                    is_subset, is_superset, is_disjoint]
                    (see: https://portcullis.readthedocs.io/en/latest/junctools.html#set)

Output:
  -- tuple of metadata and path to the merged output file, emitted as: `merged_junctions`
  -- _version.yaml containing the version string, emitted as: `version`   

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config
*/

nextflow.enable.dsl = 2

process JUNCTOOLS_SET {
    tag "${meta.id}_$mode"

    label 'process_1cpu_32gb_2h'

    input:
    tuple val(meta), path(input_files)
    val mode
  
    output:
    tuple val(meta), path('[!_]*'), emit: merged_junctions
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    junctools set \\
      $args \\
      --output ${meta.id}.merged.junctions.tab \\
      $mode \\
      $input_files

    cat <<-END_VERSION >_version.yaml
    "${task.process}":
        junctools: \$(junctools --version 2>&1 | sed -nre 's/^[^0-9]*(([0-9]+\\.)*[0-9]+).*/\\1/p')
    END_VERSION
    """
}
