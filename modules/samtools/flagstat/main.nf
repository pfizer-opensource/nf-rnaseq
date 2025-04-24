/*
Nextflow module to run samtools flagstat on an input BAM/SAM/CRAM file.

Requires: none

Params: none

Input:
  -- tuple of sample metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name)
     and a path to a bam file

Output:
  -- path to the output file, emitted as `flagstat`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process SAMTOOLS_FLAGSTAT {
    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_0'

    tag "${meta.id}"

    label 'process_6cpu_12gb_6h'

    input:
    tuple val(meta), path(bam_sorted)

    output:
    path "*.flagstat", emit: flagstat
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    """
    samtools flagstat -@ ${task.cpus} $args ${bam_sorted} > ${meta.id}.flagstat

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^samtools \\+//g; s/Using.*\$//; s/Copyright.*//')
    END_VERSION
    """
}
