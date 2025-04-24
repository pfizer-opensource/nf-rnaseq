/* Module to rebuild a STAR genome by adding novel junctions.
This is accomplished by running an alignment with a small, essentially fake,
fastq file with the novel junctions and the --sjdbInsertSave All flag. It is
essentially a trick to avoid rebuilding the entire genome with STAR genomeGenerate
and instead insert the novel junctions on-the-fly, saving the resulting genome so
that it can be used directly with any number of subsequent samples in the run. It
prevents the on-the-fly insertion from happening n times for n samples.

Requires: rebuild_genome.fastq (any small fastq file, even a single read)

Params: none

Inputs:
    -- path to the STAR SJ.out.tab file(s)
    -- path to the genome index
    -- path the GTF file

Outputs:
    -- the rebuilt genome, emitted as `genome`
    -- STAR version information, emitted as `version`

Customization:
    -- ext.args in modules.config
    -- publishDir in modules.config
*/

process STAR_REBUILD {
    container 'quay.io/biocontainers/star:2.7.1a--0'

    tag "rebuild_genome"

    label 'process_16cpu_64gb_6h'

    input:
    path junctions
    path index
    path gtf

    output:
    path './_STARgenome', emit: genome
    path '_version.yaml', emit: version

    // TODO: remove any user provided arguments that conflict with the hard-coded ones below
    // and log a warning
    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    STAR \\
      --runThreadN ${task.cpus} \\
      --genomeDir $index \\
      --sjdbGTFfile $gtf \\
      --sjdbFileChrStartEnd $junctions \\
      --readFilesIn ${projectDir}/assets/rebuild_genome.fastq \\
      --sjdbInsertSave All \\
      --outSAMtype None \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSION
    """
}
