/*
Nextflow module to run fastqc

Requires: none

Params: none

Input:
  -- tuple of metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name)
     and a path to the reads files

Output:
  -- tuple of metadata and path to a directory containing the fastqc output, emitted as `fastqc_dirs`
     the output directory will be named according to the value of `id` in the input metadata
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`
     
Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

Note: Ideally this module would request 1 thread for single-end samples and
      2 threads for paired-end. Unfortunately, dynamic labels are not technologically
      possible yet with Nextflow (https://github.com/nextflow-io/nextflow/issues/894).
      Therefore, this module simple uses the number of threads specified in the config
      files for 'process_light'. Recommend setting this to 2 for speed if resources are not
      limiting.
*/

nextflow.enable.dsl = 2

process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    tag "${meta.id}"

    label 'process_1cpu_8gb_2h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}"), emit: fastqc_dirs
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    mkdir -p ${meta.id}
    fastqc \\
        --threads ${task.cpus} \\
        --outdir ${meta.id} \\
        $args \\
        $reads

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        fastqc: \$(fastqc --version | awk '{print \$2}' | sed 's/^v//g')
    END_VERSION
    """
}
