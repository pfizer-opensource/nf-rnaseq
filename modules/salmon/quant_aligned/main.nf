/*
Nextflow module to run Salmon quantification in alignment-based mode.

Requires: none

Params: none

Input:
  -- tuple of sample metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name)
     and a path to a bam file aligned to the transcriptome (e.g., output from the STAR_ALIGN module)
  -- path to the fasta file of the transcriptome used for alignment

Output:
  -- tuple of sample metadata and the path to the Salmon output directories, named by meta.id, emitted as `salmon_dirs`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process SALMON_QUANT_ALIGNED {
    container 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'

    tag "${meta.id}"

    label 'process_6cpu_12gb_6h'

    input:
    tuple val(meta), path(bam_transcriptome)
    path transcriptome_fasta

    output:
    tuple val(meta), path("${meta.id}"), emit: salmon_dirs
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    def lib_type = args.contains('--libType') ? '' : '--libType A'
    """
    salmon quant \\
      --threads ${task.cpus} \\
      $lib_type \\
      --targets ${transcriptome_fasta} \\
      --alignments ${bam_transcriptome} \\
      --output ${meta.id} \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon \\+//g')
    END_VERSION
    """
}
