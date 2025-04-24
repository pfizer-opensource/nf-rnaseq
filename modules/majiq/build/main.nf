/*
Nextflow module to run majiq build step on STAR output bam and sj.out files.

Requires: none

Params: none

Input:
  -- majiq_settings ini file mapping samples to groups
  -- path to STAR PASS2 bams
  -- path to SAMTOOLS indices
  -- tuple of metadata and path to multipe STAR sj.out files
  -- gff3 (MAJIQ-formatted) of reference transcripts

Output:
  -- tuple of metadata and path to the output files, emitted as `majiq_out`
  -- path to splicegraph, emitted as `splicegraph`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process MAJIQ_BUILD {
    label 'process_16cpu_64gb_6h'

    input:
    path(majiq_settings)
    path(bam_sorted)
    path(bam_index)
    path(gff3)

    output:
    path("*.majiq"), emit: majiq_out
    path("*.sj"), emit: majiq_sj
    path("splicegraph.sql"), emit: splicegraph
    path("majiq.log"), emit: log
    path("_version.yaml"), emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
   
    """
    majiq build $gff3 --conf ${majiq_settings}  --nproc ${task.cpus} --output . $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        majiq: \$(echo \$(majiq -v))
    END_VERSION
    """
}

