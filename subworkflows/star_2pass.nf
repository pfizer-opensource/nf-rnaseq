/*
Run STAR_ALIGN in 2-pass mode with Portcullis filtering of pass 1 junctions.
*/

params.merge_bams = false

include { STAR_ALIGN as PASS1 } from '../modules/star/align/main' addParams(high_depth: false)
include { STAR_ALIGN as PASS2 } from '../modules/star/align/main'
include { STAR_REBUILD } from '../modules/star/rebuild/main'
include { STAR_GENOMEGENERATE } from '../modules/star/genomegenerate/main'
include { SAMTOOLS_INDEX } from '../modules/samtools/index/main'
include { SAMTOOLS_FAIDX } from '../modules/samtools/faidx/main'
include { SAMTOOLS_MERGE } from '../modules/samtools/merge/main'
include { PORTCULLIS } from '../modules/portcullis/main'
include { ANNOFILTER_JUNCTIONS } from '../modules/custom/annofilter_junctions/main'

workflow STAR_2PASS {
    take:
    reads
    index
    gtf
    junctions
    fasta

    main:

    ch_versions = Channel.empty()

    PASS1( reads, index, gtf, junctions )
    ch_versions = ch_versions.mix( PASS1.out.version.first() )
    
    SAMTOOLS_FAIDX( fasta )
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.version.first() )

    if(params.merge_bams) {
        PASS1.out.bam_sorted
            .map{ it ->
                new_meta = it[0].clone()
                new_meta.id = "${new_meta.orientation}_${new_meta.stranded}"
                tuple(new_meta, it[1])
            }
            .groupTuple()
            .set{ ch_bam_merge }
        SAMTOOLS_MERGE( ch_bam_merge )
        ch_versions = ch_versions.mix( SAMTOOLS_MERGE.out.version.first() )

        SAMTOOLS_INDEX( SAMTOOLS_MERGE.out.merged_alignment )
        ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.version.first() )
    
        PORTCULLIS( SAMTOOLS_INDEX.out.indexed_alignment, SAMTOOLS_FAIDX.out.indexed_fasta )
        ch_versions = ch_versions.mix( PORTCULLIS.out.version.first() )
        
        ANNOFILTER_JUNCTIONS( PASS1.out.junctions, PORTCULLIS.out.junctions_pass_bed.map{ it[1] }.collect(), gtf )
        ch_versions = ch_versions.mix( ANNOFILTER_JUNCTIONS.out.version.first() )
    }
    else {
        SAMTOOLS_INDEX( PASS1.out.bam_sorted )
        ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.version.first() )

        PORTCULLIS( SAMTOOLS_INDEX.out.indexed_alignment, SAMTOOLS_FAIDX.out.indexed_fasta )
        ch_versions = ch_versions.mix( PORTCULLIS.out.version.first() )

        // join channels so that STAR and portcullis files cannot become mismatched across samples
        PASS1.out.junctions
            .join(PORTCULLIS.out.junctions_pass_bed)
            .multiMap { it ->
                ch_star_junctions: tuple(it[0], it[1])
                ch_portcullis: it[2]
            }
            .set{ result }
        ANNOFILTER_JUNCTIONS( result.ch_star_junctions, result.ch_portcullis, gtf )
        ch_versions = ch_versions.mix( ANNOFILTER_JUNCTIONS.out.version.first() )
    }

    if(params.rebuild_genome) {
        STAR_REBUILD( ANNOFILTER_JUNCTIONS.out.filtered_star_junctions.map{ it -> it[1] }.collect(), index, gtf )
        ch_versions = ch_versions.mix( STAR_REBUILD.out.version.first() )
        ch_pass2_index = STAR_REBUILD.out.genome.collect()
        ch_pass2_junctions = file('NO_FILE')
    } else {
        ch_pass2_index = index
        ch_pass2_junctions = ANNOFILTER_JUNCTIONS.out.filtered_star_junctions.map{ it -> it[1] }.collect()
    }

    PASS2( reads, ch_pass2_index, gtf, ch_pass2_junctions )
    ch_versions = ch_versions.mix( PASS2.out.version.first() )

    emit:
    bam = PASS2.out.bam.ifEmpty(null)
    bam_sorted = PASS2.out.bam_sorted.ifEmpty(null)
    bam_transcriptome = PASS2.out.bam_transcriptome.ifEmpty(null)
    junctions = PASS2.out.junctions
    log_final = PASS1.out.log_final.mix(PASS2.out.log_final)
    unmapped = PASS2.out.unmapped.ifEmpty(null)
    bam_pass1 = PASS1.out.bam.ifEmpty(null)
    bam_sorted_pass1 = PASS1.out.bam_sorted.ifEmpty(null)
    bam_transcriptome_pass1 = PASS1.out.bam_transcriptome.ifEmpty(null)
    annotated_junctions = ANNOFILTER_JUNCTIONS.out.annotated_junctions
    filtered_star_junctions = ANNOFILTER_JUNCTIONS.out.filtered_star_junctions
    annofilter_log = ANNOFILTER_JUNCTIONS.out.log
    versions = ch_versions
}
