/*
Custom Nextflow module to generate an experiment.json file

The experiment.json file summarizes key metadata about the experiment and is used
to register it with REDA.

Input:
  -- reads: Channel containing a list of all read files as tuples with sample metdata (ie., (meta, reads))

Output:
  -- experiment.json, emitted as `json`

*/

import groovy.json.JsonOutput

process MAKE_EXPERIMENT_JSON {
    // label 'process_1cpu_8gb_2h'
    executor 'local'

    input:
    val reads
    tuple val(meta_genome), path(genome_fasta)

    output:
    path "experiment.json", emit: json

    exec:
    def samples = []
    def is_paired_end = []
    def stranded = []

    reads.each{
        next_sample = [
            name: it[0].id,
            order: null,
            label: null,
            meta: null,
            fastq: it[1]*.toString()
        ]
        samples.add(next_sample)
        is_paired_end.add(it[0].is_paired_end)
        stranded.add(it[0].stranded)
    }

    def experiment = [
        name: null,
        type: 'RNAseq',
        path: params.outdir,
        date: new Date().toString(),
        info: null,
        meta: [
            [
                name: 'runName',
                value: workflow.runName,
                label: 'Nextflow Run Name'
            ],
            [
                name: 'sessionId',
                value: workflow.sessionId,
                label: 'Nextflow Session ID'
            ]
        ],
        contacts: [
            [
            name: null,
            role: 'analyst',
            login: workflow.userName,
            email: null,
            group: null
            ]
        ],
        project: null,
        organism: meta_genome.species,
        stranded: stranded.unique().join(','),
        is_paired_end: null,
        is_paired_end: is_paired_end.unique().join(','),
        samples: samples
    ]
    def json = JsonOutput.toJson(experiment)
    outfile = task.workDir.resolve('experiment.json')
    outfile.text = json
}
