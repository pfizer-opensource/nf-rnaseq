/*
Nextflow module to run Portcullis. This is essentially `portcullis full` but implemented step-by-step.
This allows finer control over optional arguments for each of the steps.

Requires: samtools

Params:
  -- portcullis_maxthreads: (integer, default = 8) maximum number of threads to use for Portcullis;
                            Portcullis may use more memory with little gain in run-time past about 8 threads

Input:
  -- tuple of metadata (a Groovy map containing, at minimum, keys called `id` (a sample identifier), `stranded` 
     (strandedness of the library prep), and `orientation` (orientation of the reads, e.g., FR) and a path to BAM file(s)
     NOTE: it is likely more efficient to manually merge multiple bam files (e.g., SAMTOOLS_MERGE) first using appropriate
     resource settings rather than allowing Portcullis to call samtools merge with the typical 4 or 8 threads that are optimal
     for Portcullis
  -- fasta representing the genome used to align the BAM file

Output (all are a tuple of sample metadata and the following paths):
  -- path to the portcullis junc output in tab format, emitted as: `junctions_tab`
  -- path to the portcullis junc output in bed format, emitted as: `junctions_bed`
  -- path to the portcullis filt passing junctions in tab format (post-MT/scaffold
     filtering if enabled), emitted as: `junctions_pass_tab`
  -- path to the portcullis filt passing junctions output in bed format, emitted as: `junctions_pass_bed`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`
  -- optional outputs (if enabled by optional arguments):
    -- path to the portcullis filt failing junctions in tab format, emitted as: `junctions_fail_tab`
    -- path to the portcullis filt failing junctions in bed format, emitted as: `junctions_fail_bed`      
    -- path to portcullis filt passing junctions as exon-based gff, emitted as: `exon_pass_gff`
    -- path to portcullis filt failing junctions as exon-based gff, emitted as: `exon_fail_gff`      
    -- path to portcullis filt passing junctions as intron-based gff, emitted as: `intron_pass_gff`
    -- path to portcullis filt failing junctions as intron-based gff, emitted as: `intron_fail_gff`      

Customization:
    IMPORTANT: Changing -o/--output in module.config will likely break this module! Strongly advise to leave defaults!
    ext.args.prep in module.config to add optional arguments to portcullis prep
    ext.args.junc in module.config to add optional arguments to portcullis junc
    ext.args.filt in module.config to add optional arguments to portcullis filt
    publishDir in module.config
*/

nextflow.enable.dsl = 2

params.portcullis_maxthreads = 8

process PORTCULLIS {
    tag "${meta.id}"

    label 'process_4cpu_256gb_6h'

    input:
    tuple val(meta), path(bam), path(bam_index)
    tuple path(fasta), path(fasta_index)

    output:
    tuple val(meta), path('**junc/*.junctions.tab'), emit: junctions_tab
    tuple val(meta), path('**junc/*.junctions.bed'), emit: junctions_bed
    tuple val(meta), path('**.pass.junctions.tab'), emit: junctions_pass_tab
    tuple val(meta), path('**.pass.junctions.bed'), emit: junctions_pass_bed
    tuple val(meta), path('**.fail.junctions.tab'), emit: junctions_fail_tab, optional: true
    tuple val(meta), path('**.fail.junctions.bed'), emit: junctions_fail_bed, optional: true
    tuple val(meta), path('**.pass.junctions.exon.gff'), emit: exon_pass_gff, optional: true
    tuple val(meta), path('**.fail.junctions.exon.gff'), emit: exon_fail_gff, optional: true
    tuple val(meta), path('**.pass.junctions.intron.gff'), emit: intron_pass_gff, optional: true
    tuple val(meta), path('**.fail.junctions.intron.gff'), emit: intron_fail_gff, optional: true
    path "_version.yaml", emit: version

    script:
    def prep_args = task.ext.args.prep ? task.ext.args.prep.join(' \\\n      ').trim() : ''
    def junc_args = task.ext.args.junc ? task.ext.args.junc.join(' \\\n      ').trim() : ''
    def filt_args = task.ext.args.filt ? task.ext.args.filt.join(' \\\n      ').trim() : ''

    def num_threads = Math.max(task.cpus, params.portcullis_maxthreads)

    // other tools use strandedness information using a different vocabulary so a translator is necessary
    // check if strandedness has been defined using Portcullis terminology first
    // if not, attempt to convert using a commonly used vocabulary
    // if that fails, Portcullis will emit an error
    def valid_strandedness = ["firststrand", "secondstrand", "unstranded", "UNKNOWN"]
    if(valid_strandedness.contains(meta.stranded)) {
      strand_arg = meta.stranded
    }
    else {
        switch(meta.stranded) {
            case "forward":
            strand_arg = "firststrand"
            break
            case "reverse":
            strand_arg = "secondstrand"
            break
            case "unstranded":
            strand_arg = "unstranded"
            break
            case "none":
            strand_arg = "unstranded"
            break
            default:
            log.warn("strandedness must be one of ['forward', 'reverse', 'none', 'unstranded'], not ${meta.stranded}: using --stranded UNKNOWN")
            strand = "UNKNOWN"
          }
        }

    """
    portcullis prep \\
      --threads $num_threads \\
      --output ${meta.id}/prep \\
      $prep_args \\
      $fasta \\
      $bam

    portcullis junc \\
      --threads $num_threads \\
      --output ${meta.id}/junc/${meta.id}.portcullis \\
      $junc_args \\
      --strandedness $strand_arg \\
      --orientation ${meta.orientation} \\
      ${meta.id}/prep

    portcullis filter \\
      --threads $num_threads \\
      --output ${meta.id}/${meta.id}.portcullis \\
      $filt_args \\
      ${meta.id}/prep \\
      ${meta.id}/junc/${meta.id}.portcullis.junctions.tab

    cat <<-END_VERSION >_version.yaml
    "${task.process}":
        portcullis: \$(portcullis --version 2>&1 | sed -nre 's/^[^0-9]*(([0-9]+\\.)*[0-9]+).*/\\1/p')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^samtools \\+//g; s/Using.*\$//; s/Copyright.*//')
        bcftools: \$(bcftools --version 2>&1 | grep bcftools | sed -nre 's/^[^0-9]*(([0-9]+\\.)*[0-9]+).*/\\1/p')
        htslib: \$(bcftools --version 2>&1 | grep htslib | sed -nre 's/^[^0-9]*(([0-9]+\\.)*[0-9]+).*/\\1/p')
    END_VERSION
    """
}