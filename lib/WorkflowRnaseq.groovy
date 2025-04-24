//
// This file holds several functions specific to the workflow/rnaseq.nf
//

class WorkflowRnaseq {

    //
    // Check and validate parameters
    //
    public static void initialize(params, log) {
        // Verify that provided input directory provides valid pipeline input
        validateInput(params, log)

        // Verify that provided outdir directory is writeable
        validateOutdirWritable(params, log)

        // Verify that provided genome_registry path is a valid registry
        genomeRegistryExistsError(params, log)

        // TODO: add validation for optional STAR junctions file
        // schema: chr, start, end, strand

        // Check which RSeQC modules we are running
        // def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
        // if ((valid_params['rseqc_modules'] + rseqc_modules).unique().size() != valid_params['rseqc_modules'].size()) {
        //     log.error "Invalid option: ${params.rseqc_modules}. Valid options for '--rseqc_modules': ${valid_params['rseqc_modules'].join(', ')}"
        //     System.exit(1)
        // }
    }

    private static void validateInput(params, log) {
        def input_directory = new File(params.input)
        if (!input_directory.isDirectory()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Input directory not found at '${params.input}'.\n" +
                "  Check path to input directory.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!input_directory.canRead()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Input directory at '${params.input}' cannot be read.\n" +
                "  Check permissions on input directory.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        // verify fastq directory availabl and readable
        def fastq_directory = new File(params.fastq_dir)
        if (!fastq_directory.isDirectory()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Input directory does not contain a subdirectory named 'fastq'.\n" +
                "  Check that fastq files are provided in a folder named 'fastq' in the input directory.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!fastq_directory.canRead()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Fastq directory at '${fastq_directory}' cannot be read.\n" +
                "  Check permissions on fastq directory.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        // verify sample sheet is available and readable
        def samplesheet = new File(params.input, 'SampleSheet.csv')
        if (!samplesheet.isFile()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  SampleSheet.csv not found in '${params.input}'.\n" +
                "  Check that input directory contains a sample sheet named 'SampleSheet.csv'.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!samplesheet.canRead()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Sample sheet at '${samplesheet}' cannot be read.\n" +
                "  Check permissions on SampleSheet.csv.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        // if MAJIQ is being run, verify required input file
        if (params.run_majiq) {
            def majiq_contrasts_csv = new File(params.input, 'majiq_contrasts.csv')
            if (!majiq_contrasts_csv.isFile()) {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  params.run_majiq = true but required input file majiq_contrasts.csv not found.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }
          
            if (!majiq_contrasts_csv.canRead()) {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  params.run_majiq = true but required input file majiq_contrasts.csv is not readable.\n" +
                    "  Check that majiq_contrasts.csv is valid with proper read permissions.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }
        }
    }

    private static void validateOutdirWritable(params, log) {
        def output_directory = new File(params.outdir)
        if (output_directory.parent != params.input) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Output directory '${params.outdir}' is in different location than input directory '${params.input}.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }

        try {
            output_directory.mkdir()
        } catch(Exception ex) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Output directory cannot be created at '${params.outdir}'.\n" +
                "  Check appropriate permissions to write to output path.\n" +
                "  ${ex.toString()}" +
                "  ${ex.getMessage()}" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    private static void genomeRegistryExistsError(params, log) {
        // def registry_paths = [params.genome_registry, 'genomes']
        // def genome_registry = new File(registry_paths.join(File.separator))
        def genome_registry = new File(params.genome_registry, 'genomes')
        if (!genome_registry.exists()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome registry not found at '${params.genome_registry}'.\n" +
                "  Check path to genome registry.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    // //
    // // Function to check whether biotype field exists in GTF file
    // //
    // public static Boolean biotypeInGtf(gtf_file, biotype, log) {
    //     def hits = 0
    //     gtf_file.eachLine { line ->
    //         def attributes = line.split('\t')[-1].split()
    //         if (attributes.contains(biotype)) {
    //             hits += 1
    //         }
    //     }
    //     if (hits) {
    //         return true
    //     } else {
    //         log.warn "=============================================================================\n" +
    //             "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
    //             "  Biotype QC will be skipped to circumvent the issue below:\n" +
    //             "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
    //             "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
    //             "==================================================================================="
    //         return false
    //     }
    // }

    // //
    // // Function to generate an error if contigs in genome fasta file > 512 Mbp
    // //
    // public static void checkMaxContigSize(fai_file, log) {
    //     def max_size = 512000000
    //     fai_file.eachLine { line ->
    //         def lspl  = line.split('\t')
    //         def chrom = lspl[0]
    //         def size  = lspl[1]
    //         if (size.toInteger() > max_size) {
    //             log.error "=============================================================================\n" +
    //                 "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
    //                 "  ${chrom}: ${size}\n\n" +
    //                 "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
    //                 "  Please see:\n" +
    //                 "  https://github.com/nf-core/rnaseq/issues/744\n" +
    //                 "============================================================================="
    //             System.exit(1)
    //         }
    //     }
    // }

    // //
    // // Function that parses TrimGalore log output file to get total number of reads after trimming
    // //
    // public static Integer getTrimGaloreReadsAfterFiltering(log_file) {
    //     def total_reads = 0
    //     def filtered_reads = 0
    //     log_file.eachLine { line ->
    //         def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
    //         def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
    //         if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
    //         if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    //     }
    //     return total_reads - filtered_reads
    // }

    // //
    // // Function that parses and returns the alignment rate from the STAR log output
    // //
    // public static ArrayList getStarPercentMapped(params, align_log) {
    //     def percent_aligned = 0
    //     def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    //     align_log.eachLine { line ->
    //         def matcher = line =~ pattern
    //         if (matcher) {
    //             percent_aligned = matcher[0][1].toFloat()
    //         }
    //     }

    //     def pass = false
    //     if (percent_aligned >= params.min_mapped_reads.toFloat()) {
    //         pass = true
    //     }
    //     return [ percent_aligned, pass ]
    // }

    // //
    // // Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
    // //
    // public static ArrayList getInferexperimentStrandedness(inferexperiment_file, cutoff=30) {
    //     def sense        = 0
    //     def antisense    = 0
    //     def undetermined = 0
    //     inferexperiment_file.eachLine { line ->
    //         def undetermined_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
    //         def se_sense_matcher     = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
    //         def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
    //         def pe_sense_matcher     = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
    //         def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/
    //         if (undetermined_matcher) undetermined = undetermined_matcher[0][1].toFloat() * 100
    //         if (se_sense_matcher)     sense        = se_sense_matcher[0][1].toFloat() * 100
    //         if (se_antisense_matcher) antisense    = se_antisense_matcher[0][1].toFloat() * 100
    //         if (pe_sense_matcher)     sense        = pe_sense_matcher[0][1].toFloat() * 100
    //         if (pe_antisense_matcher) antisense    = pe_antisense_matcher[0][1].toFloat() * 100
    //     }
    //     def strandedness = 'unstranded'
    //     if (sense >= 100-cutoff) {
    //         strandedness = 'forward'
    //     } else if (antisense >= 100-cutoff) {
    //         strandedness = 'reverse'
    //     }
    //     return [ strandedness, sense, antisense, undetermined ]
    // }

    // //
    // // Get workflow summary for MultiQC
    // //
    // public static String paramsSummaryMultiqc(workflow, summary) {
    //     String summary_section = ''
    //     for (group in summary.keySet()) {
    //         def group_params = summary.get(group)  // This gets the parameters of that particular group
    //         if (group_params) {
    //             summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
    //             summary_section += "    <dl class=\"dl-horizontal\">\n"
    //             for (param in group_params.keySet()) {
    //                 summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
    //             }
    //             summary_section += "    </dl>\n"
    //         }
    //     }

    //     String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
    //     yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
    //     yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    //     yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    //     yaml_file_text        += "plot_type: 'html'\n"
    //     yaml_file_text        += "data: |\n"
    //     yaml_file_text        += "${summary_section}"
    //     return yaml_file_text
    // }//
    // Exit pipeline if incorrect --genome_registry key provided
    //

    // //
    // // Print a warning if using GRCh38 assembly from igenomes.config
    // //
    // private static void ncbiGenomeWarn(log) {
    //     log.warn "=============================================================================\n" +
    //         "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
    //         "  Biotype QC will be skipped to circumvent the issue below:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
    //         "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
    //         "==================================================================================="
    // }

    // //
    // // Print a warning if using a UCSC assembly from igenomes.config
    // //
    // private static void ucscGenomeWarn(log) {
    //     log.warn "=============================================================================\n" +
    //         "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
    //         "  Biotype QC will be skipped to circumvent the issue below:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
    //         "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
    //         "==================================================================================="
    // }

    // //
    // // Print a warning if both GTF and GFF have been provided
    // //
    // private static void gtfGffWarn(log) {
    //     log.warn "=============================================================================\n" +
    //         "  Both '--gtf' and '--gff' parameters have been provided.\n" +
    //         "  Using GTF file as priority.\n" +
    //         "==================================================================================="
    // }

    // //
    // // Print a warning if using '--transcript_fasta'
    // //
    // private static void transcriptsFastaWarn(log) {
    //     log.warn "=============================================================================\n" +
    //         "  '--transcript_fasta' parameter has been provided.\n" +
    //         "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" +
    //         "  Please see:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/753\n" +
    //         "==================================================================================="
    // }

    // //
    // // Print a warning if --skip_alignment has been provided
    // //
    // private static void skipAlignmentWarn(log) {
    //     log.warn "=============================================================================\n" +
    //         "  '--skip_alignment' parameter has been provided.\n" +
    //         "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" +
    //         "==================================================================================="
    // }

    // //
    // // Print a warning if using '--aligner star_rsem' and '--with_umi'
    // //
    // private static void rsemUmiError(log) {
    //     log.error "=============================================================================\n" +
    //         "  When using '--aligner star_rsem', STAR is run by RSEM itself and so it is\n" +
    //         "  not possible to remove UMIs before the quantification.\n\n" +
    //         "  If you would like to remove UMI barcodes using the '--with_umi' option\n" +
    //         "  please use either '--aligner star_salmon' or '--aligner hisat2'.\n" +
    //         "============================================================================="
    //     System.exit(1)
    // }

    // //
    // // Print a warning if using '--aligner star_rsem' and providing both '--rsem_index' and '--star_index'
    // //
    // private static void rsemStarIndexWarn(log) {
    //     log.warn "=============================================================================\n" +
    //         "  When using '--aligner star_rsem', both the STAR and RSEM indices should\n" +
    //         "  be present in the path specified by '--rsem_index'.\n\n" +
    //         "  This warning has been generated because you have provided both\n" +
    //         "  '--rsem_index' and '--star_index'. The pipeline will ignore the latter.\n\n" +
    //         "  Please see:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/568\n" +
    //         "==================================================================================="
    // }

    // //
    // // Print a warning if using '--additional_fasta' and '--<ALIGNER>_index'
    // //
    // private static void additionaFastaIndexWarn(index, log) {
    //     log.warn "=============================================================================\n" +
    //         "  When using '--additional_fasta <FASTA_FILE>' the aligner index will not\n" +
    //         "  be re-built with the transgenes incorporated by default since you have \n" +
    //         "  already provided an index via '--${index}_index <INDEX>'.\n\n" +
    //         "  Set '--additional_fasta <FASTA_FILE> --${index}_index false --save_reference' to\n" +
    //         "  re-build the index with transgenes included and the index will be saved in\n" +
    //         "  'results/genome/index/${index}/' for re-use with '--${index}_index'.\n\n" +
    //         "  Ignore this warning if you know that the index already contains transgenes.\n\n" +
    //         "  Please see:\n" +
    //         "  https://github.com/nf-core/rnaseq/issues/556\n" +
    //         "==================================================================================="
    // }
}