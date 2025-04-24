/*
Module to run STAR in alignment mode

Requires: none

Params:
  -- high_depth: (boolean) set to true to switch process label to process_20cpu_48gb_12h providing a separate resource request
                 for deeply sequenced runs
  -- dev: (boolean) run in dev mode; used to specify smaller executor resource requirements when
          testing with small test data sets
Inputs:
  -- (meta, reads): tuple of metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name)
     and a path to the reads files
  -- index: path to the genome index (e.g., as created by STAR_GENOMEGENERATE)
  -- gtf: path to the GTF file
  -- junctions: optional path to all junction files to include (e.g., 2nd pass alignment)
     NOTE: if no junctions to include, must pass a file named 'NO_FILE' (workaround for optional input
     suggested at nextflow-io.github.io/patterns/optional-input)

Output:
  NOTE: many outputs are optional because the specific run arguments will determine which ones are
        produced (i.e., optional arguments supplied in module.config)
  -- tuple of metadata and unsorted bam file, emitted as `bam`
  -- tuple of metadata and sorted bam file, emitted as `bam_sorted`
  -- tuple of metadata and bam of transcriptome alignment, emitted as `bam_transcriptome`
  -- tuple of metadata and the junctions file, emitted as `junctions`
  -- tuple of metadata and the unmapped reads file(s), emitted as `unmapped`
  -- tuple of metadata and the Log.final.out, emitted as `log_final`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in modules.config to add optional arguments
    ext.suffix in modules.config to add a custom suffix to meta.id; this is used as the --outFileNamePrefix argument
    publishDir in modules.config
*/

process STAR_ALIGN {
    container 'quay.io/biocontainers/star:2.7.1a--0'

    tag "${meta.id}"

    label params.dev ? 'process_8cpu_32gb_2h' : params.high_depth ? 'process_20cpu_48gb_12h' : 'process_16cpu_64gb_6h'

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    path junctions

    output:
    tuple val(meta), path('*.Aligned.out.bam'), optional: true, emit: bam
    tuple val(meta), path('*.Aligned.sortedByCoord.out.bam'), optional: true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional: true, emit: bam_transcriptome
    tuple val(meta), path('*.SJ.out.tab'), emit: junctions
    tuple val(meta), path('*.Unmapped.out.mate{1,2}.fastq.gz'), optional: true, emit: unmapped
    tuple val(meta), path('*.Log.final.out'), emit: log_final
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    def prefix = task.ext.suffix ? "${meta.id}.${task.ext.suffix}" : "${meta.id}"
    def read_files_command = args.contains('--readFilesCommand') ? '' : '--readFilesCommand "gzip -dc"'
    def out_sam_type = args.contains('--outSAMtype') ? '' : '--outSAMtype BAM Unsorted'
    def sjdbfile = junctions.name != 'NO_FILE' ? "--sjdbFileChrStartEnd $junctions" : ''
    def gzip_reads = args.contains('--outReadsUnmapped Fastx') ? 'gzip *Unmapped* && for file in *Unmapped*; do mv "$file" "${file:0: -3}.fastq.gz"; done' : ''
    """
    STAR \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${index} \\
      --readFilesIn ${reads} \\
      --outFileNamePrefix $prefix. \\
      --sjdbGTFfile ${gtf} \\
      $read_files_command \\
      $sjdbfile \\
      $out_sam_type \\
      $args

    $gzip_reads

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSION
    """
}
