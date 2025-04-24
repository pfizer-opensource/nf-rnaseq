/*
Nextflow module to run samtools merge on a BAM/SAM/CRAM file.

Requires: none

Params: none

Input:
  -- tuple of metadata and path to multiple BAM/SAM/CRAM files

Output:
  -- tuple of metadata and path to the merged output file, emitted as `merged_alignment`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process SAMTOOLS_MERGE {
    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_0'

    tag "${meta.id}"

    label 'process_16cpu_64gb_6h'

    input:
    tuple val(meta), path(sorted_alignments)

    output:
    tuple val(meta), path("[!_]*"), emit: merged_alignment
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    def outname = args.contains('-o') ? '' : "${meta.id}.merged.bam"
    """
    samtools merge -@ ${task.cpus} $args $outname ${sorted_alignments}

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^samtools \\+//g; s/Using.*\$//; s/Copyright.*//')
    END_VERSION
    """
}
