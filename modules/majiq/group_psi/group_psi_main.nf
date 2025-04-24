/*
Nextflow module to run majiq psi step on groups of samples

Requires: none

Params: none

Input:
  -- path to MAJIQ_BUILD .majiq files

Output:
  -- path to the per-group psi output files, emitted as `group_psi`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process MAJIQ_GROUP_PSI {
    label 'process_1cpu_8gb_2h'

    input:
    val(majiq_group)
    path(majiq_out)

    output:
    path("*.psi.tsv"), emit: group_psi_tsv
    path("*.psi.voila"), emit: group_psi_voila
    path("*psi_majiq.log"), emit: log
    path("_version.yaml"), emit: version

    script:
    def args = task.ext.args ?: ''
    """   
    majiq psi ${majiq_out} --nproc ${task.cpus} --output . --name ${majiq_group}
    mv psi_majiq.log ${majiq_group}.psi_majiq.log
    
    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        majiq: \$(echo \$(majiq -v))
    END_VERSION
    """
}

