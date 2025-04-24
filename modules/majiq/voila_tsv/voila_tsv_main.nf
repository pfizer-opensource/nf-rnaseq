/*
Nextflow module to run majiq build step on STAR output bam and sj.out files.

Requires: none

Params: none

Input:
  -- path to majiq splicegraph.sql from majiq build step
  -- path to deltapsi.voila 

Output:
  -- path to tsv file for voila view
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process VOILA_TSV {
    label 'process_1cpu_8gb_2h'

    input:
    path(splicegraph)
    path(deltapsi)

    output:
    path("*.deltapsi.voila.tsv"), emit: deltapsi_voila_tsv
    path("_version.yaml"), emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    
    voila tsv ${splicegraph} ${deltapsi} --file-name ${deltapsi}.tsv $args
    
    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        majiq_voila: \$(echo \$(voila -v))
    END_VERSION
    """
}

