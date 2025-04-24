/*
Nextflow module to run Picard CollectRnaSeqMetrics

Requires: none

Params:
  -- dev: (boolean) run in dev mode; used to specify smaller executor resource requirements when
    testing with small test data sets

Input:
  -- tuple of sample metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name)
     and a path to the bam file
  -- a refflat file corresponding to the exact genome used in alignment that generated the input bam file
     see: https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard-#--REF_FLAT
     NOTE: this can be generated using the gtf_to_refflat module
  -- a file containing the rRNA intervals in the exact genome used in the alignment that generated the input bam file
     see: http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html
     NOTE: this can be generated using the make_ribosomal_intervals module

Output:
  -- tuple of sample metadata and path to the output report file
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

nextflow.enable.dsl = 2

process PICARD_COLLECTRNASEQMETRICS {
    container 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'

    tag "${meta.id}"

    label params.dev ? 'process_8cpu_32gb_2h' : 'process_6cpu_12gb_6h'

    input:
    tuple val(meta), path(bam)
    path refflat
    path ribosomal_intervals

    output:
    tuple val(meta), path("${meta.id}.RNA_Metrics"), emit: picard_rna_metrics
    path "_version.yaml", emit: version

    script:

    String strand
    switch(meta.stranded) {
        case 'forward':
            strand = 'FIRST_READ_TRANSCRIPTION_STRAND'
            break
        case 'reverse':
            strand = 'SECOND_READ_TRANSCRIPTION_STRAND'
            break
        case 'unstranded':
            strand = 'NONE'
            break
        case 'none':
            strand = 'NONE'
            break
        default:
            strand = 'NONE'
    }
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    picard \\
      -Xms4096m \\
      -Xmx4096m \\
      -XX:ParallelGCThreads=${task.cpus} \\
      CollectRnaSeqMetrics \\
      --INPUT $bam \\
      --OUTPUT "${meta.id}.RNA_Metrics" \\
      --STRAND_SPECIFICITY $strand \\
      --REF_FLAT $refflat \\
      --RIBOSOMAL_INTERVALS $ribosomal_intervals \\
      --VERBOSITY DEBUG \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        picard: \$(picard CollectRnaSeqMetrics --version 2>&1 | tail -1 | grep -Eo '[0-9]+(\\.[0-9]+)*')
    END_VERSION
    """
}
