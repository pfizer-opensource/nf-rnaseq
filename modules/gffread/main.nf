/*
Nextflow module to run gffread

Requires: none

Params: none

Input:
  -- gtf: path to a GTF file
  -- optional_fasta: optional path to a gzipped fasta file (default of 'NO_FILE' indicates no fasta file)

Output:
  -- optional path to a GTF/GFF file, emitted as `gtf`
  -- optional path to a fasta file, emitted as `fasta`
  -- optional path to a bed file, emitted as `bed`
  -- optional path to a tsv file, emitted as `tsv`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process GFFREAD {
    tag "$gtf"

    label 'process_1cpu_8gb_2h'

    input:
    path gtf
    path optional_fasta

    output:
    path '*.{gtf,gff}', optional: true, emit: gtf
    path '*.{fa,fasta}', optional: true, emit: fasta
    path '*.bed', optional: true, emit: bed
    path '*.tsv', optional: true, emit: tsv
    path "_version.yaml", emit: version

    script:
    def fasta = optional_fasta.name != 'NO_FILE' ? "-g $optional_fasta" : ''
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim(): ''
    """
    gffread \\
        $fasta \\
        $args \\
        <(zcat $gtf)

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gffread: \$(gffread --version)
    END_VERSION
    """
}
