/*
Nextflow DSL2 module to compress files using pigz.
*/

process PIGZ {
    tag "$uncompressed"

    label 'process_8cpu_24gb_30min'

    stageInMode 'copy'

    input:
    tuple val(meta), path(uncompressed)

    output:
    tuple val(meta), path("*.gz"), emit: compressed
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    """
    pigz -k $args $uncompressed

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1 | sed "s/^pigz //"))
    END_VERSION
    """
}
