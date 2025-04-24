/*
Run STAR_ALIGN in 2-pass mode.
*/

include { STAR_ALIGN as PASS1 } from '../modules/star/align/main' addParams(star_mode: 'pass1')
include { STAR_ALIGN as PASS2 } from '../modules/star/align/main' addParams(star_mode: 'pass2')

workflow STAR_2PASS {
    take:
    reads
    index
    gtf
    junctions

    main:
    PASS1( reads, index, gtf, junctions )
    PASS2( reads, index, gtf, PASS1.out.junctions.collect{ it[1] } )

    emit:
    bam = PASS2.out.bam.ifEmpty(null)
    bam_sorted = PASS2.out.bam_sorted.ifEmpty(null)
    junctions = PASS2.out.junctions
    log_final = PASS2.out.log_final
    unmapped = PASS2.out.unmapped.ifEmpty(null)
//    tuple val(meta), path('*.fastq.gz'), optional: true, emit: fastq
    bam_transcriptome = PASS2.out.bam_transcriptome.ifEmpty(null)
}
