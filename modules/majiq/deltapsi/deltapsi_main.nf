/*
Nextflow module to run majiq deltapsi step using majiq build outputs files.

Requires: none

Params: none

Input:
  -- tuple of group names
  -- paths to majiq_reference_samples
  -- paths to majiq_experiment_samples
  -- path to splicegraph from majiq build

Output:
  -- path to the per-comparison deltapsi.tsv output files, emitted as `sample_deltapsi_tsv`
  -- path to the per-comparison deltapsi.voila output files, emitted as `sample_deltapsi_voila`
  -- path to the deltapsi log file, emitted as `log`
  -- path to tsv file for voila view
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process MAJIQ_DELTAPSI {
    label 'process_8cpu_16gb_30min'

    input:
    tuple val(contrast_name), val(reference_group), val(experiment_group)
    path(majiq_reference_samples)
    path(majiq_experiment_samples)
    path(splicegraph)

    output:
    tuple val(contrast_name), path("*.deltapsi.tsv"), emit: deltapsi_tsv
    path("*.deltapsi.voila"), emit: deltapsi_voila
    path("*.deltapsi_majiq.log"), emit: log
    path("*.deltapsi.voila.tsv"), emit: deltapsi_voila_tsv
    path("_version.yaml"), emit: version

    script:
    def args1 = task.ext.args1 ? task.ext.args1.join(' \\\n      ').trim() : ''
    def args2 = task.ext.args2 ? task.ext.args2.join(' \\\n      ').trim() : ''
    """
    
    majiq deltapsi -grp1 ${majiq_reference_samples} -grp2 ${majiq_experiment_samples} --nproc ${task.cpus} --output . --names ${reference_group} ${experiment_group} $args1
    mv deltapsi_majiq.log ${reference_group}-${experiment_group}.deltapsi_majiq.log
    
    voila tsv ${splicegraph} ${reference_group}-${experiment_group}.deltapsi.voila -f ${reference_group}-${experiment_group}.deltapsi.voila.tsv $args2
    
    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        majiq: \$(echo \$(majiq -v))
    END_VERSION
    """
}

