/*
Nextflow module to run samtools index on a BAM/CRAM or BGZF-compressed SAM file.

Requires: none

Params: none

Input:
  -- tuple of metadata map and path to a BAM/CRAM or BGZF-compressed SAM file

Output:
  -- tuple of metadata, path to the input BAM/SAM/CRAM file, and path to the index emitted as `indexed_alignment`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process SAMTOOLS_INDEX {
    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_0'

    tag "$bam_sorted"

    label 'process_6cpu_12gb_6h'

    input:
    tuple val(meta), path(bam_sorted)

    output:
    tuple val(meta), path("*.{bam,cram,sam.gz}", includeInputs: true), path("*.{bai,csi,crai}"), emit: indexed_alignment
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    """
    samtools index -@ ${task.cpus} $args ${bam_sorted}

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^samtools \\+//g; s/Using.*\$//; s/Copyright.*//')
    END_VERSION
    """
}
