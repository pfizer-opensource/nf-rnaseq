/*
Nextflow module to run samtools faidx on a fasta file.
TODO: Add optional regions. Currently only generates an index.

Requires: none

Params: none

Input:
  -- path to a fasta file

Output:
  -- tuple of path to the indexed fasta file and path to the .fai output file, emitted as `fasta_indexed`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process SAMTOOLS_FAIDX {
    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_0'

    tag "$fasta"

    label 'process_1cpu_32gb_2h'

    input:
    path fasta

    output:
    tuple path("*.{fa,fasta}", includeInputs: true), path("*.fai"), emit: indexed_fasta
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    """
    samtools faidx $args ${fasta}

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^samtools \\+//g; s/Using.*\$//; s/Copyright.*//')
    END_VERSION
    """
}
