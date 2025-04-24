/*
Nextflow module to run majiq build step on STAR output bam and sj.out files.

Requires: none

Params: none

Input:
  -- path to MAJIQ_BUILD .majiq files

Output:
  -- path to the per-sample psi output files, emitted as `sample_psi`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process MAJIQ_SAMPLE_PSI {
    label 'process_1cpu_8gb_2h'

    input:
    path(majiq_out)

    output:
    path("*.psi.tsv"), emit: sample_psi_tsv
    path("*.psi.voila"), emit: sample_psi_voila
    path("*psi_majiq.log"), emit: log
    path("_version.yaml"), emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    majiq psi ${majiq_out} --nproc ${task.cpus} --output . --name tempname $args
    mv tempname.psi.tsv ${majiq_out}.psi.tsv
    mv tempname.psi.voila ${majiq_out}.psi.voila
    
    mv psi_majiq.log ${majiq_out}.psi_majiq.log
    
    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        majiq: \$(echo \$(majiq -v))
    END_VERSION
    """
}

