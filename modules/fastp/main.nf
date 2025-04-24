/*
Nextflow module to run Fastp. 
Requires: none

Params: none

Input:
  -- tuple of metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name and
     an entry called 'is_paired_end' indicating whether reads are paired vs. single), and a path to the reads files

Output (all are a tuple of sample metadata and the following paths):
  -- path to the fastq reads output, emitted as: `fastp_reads`
  -- path to fastp.html, emitted as: `fastp_html`
  -- path to fastp.json, emitted as: `fastp_json` 
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in modules.config to add optional arguments fastp
    publishDir in modules.config
*/

nextflow.enable.dsl = 2

process FASTP {
    container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'

    tag "${meta.id}"

    label 'process_6cpu_12gb_6h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fq.gz'), emit: fastp_reads
    tuple val(meta), path('*.fastp.html'), emit: fastp_html
    tuple val(meta), path('*.fastp.json'), emit: fastp_json
    path "_version.yaml", emit: version

    script:
    def in_arg = meta.is_paired_end ? "--in1 ${reads[0]} --in2 ${reads[1]}": "--in1 $reads"
    def out_arg = meta.is_paired_end ? "--out1 ${meta.id}_trimmed_R1.fq.gz --out2 ${meta.id}_trimmed_R2.fq.gz": "--out1 ${meta.id}.trimmed.fq.gz"
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''   
    """
    fastp \\
      $in_arg \\
      $out_arg \\
      --html ${meta.id}.fastp.html \\
      --json ${meta.id}.fastp.json \\
      --report_title ${meta.id} \\
      --thread ${task.cpus} \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)*')
    END_VERSION
    """
}
