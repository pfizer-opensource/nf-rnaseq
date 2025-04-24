/*
Custom module to ANNOtate Portcullis junctions and to FILTER a STAR SJ.out.tab
file using the results

Requires: annotate_junctions.py

Params: 
    -- merge_bams (boolean) combine bam files for all runs

Inputs:
    -- tuple of metadata (a Groovy map containing, at minimum, and entry called `id` containing the sample name)
       and a path to the STAR SJ.out.tab file
    -- path to the Portcullis bed output file (e.g., portcullis.pass.bed)
    -- path to the GTF file

Outputs:
    -- tuple metadata and path to the filtered STAR SJ.out.tab file, emitted as `filtered_star_junctions`
    -- tuple metadata and path to the annotated junctions file, emitted as `annotated_junctions`
    -- path to the annotate_junctions.py log file, emitted as `log`
    -- version info yaml file, emitted as `version`
*/

nextflow.enable.dsl=2

process ANNOFILTER_JUNCTIONS {
    tag "${meta.id}"

    label params.merge_bams ? 'process_1cpu_24gb_2hr' : 'process_1cpu_8gb_2h'

    input:
    tuple val(meta), path(star_junctions)
    path portcullis_bed
    path gtf

    output:
    tuple val(meta), path('*SJ.out.filtered.tab'), emit: filtered_star_junctions
    tuple val(meta), path('*annotated.tsv'), emit: annotated_junctions
    path '*.log', emit: log
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    annofilter_junctions.py \\
      ${portcullis_bed} \\
      ${gtf} \\
      --star ${star_junctions} \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        annofilter_junctions.py: \$(echo \$(annofilter_junctions.py --version) | sed 's/annofilter_junctions.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
