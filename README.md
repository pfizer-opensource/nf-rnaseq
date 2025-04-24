# Nextflow RNAseq
## Complete pipeline for short-read RNAseq processing including STAR 2-pass and alternative splicing analysis
This pipeline implements the popular STAR-Salmon workflow with the option for either 1-pass and 2-pass alignment mode. The 2-pass mode is embued with additional tools and options that aim to minimize incorporation of false junctions in the sencond pass. It is expected that novel junctions will lead some previously uniquely mapping reads in the first pass to multimap in the second pass. Provided these novel junctions are real, this is an imporovement over 1-pass alignment. However, if a uniquely mapping read is multimapped to a false junction in pass 2, the alignment has actually become worse in 2-pass mode. With that concern in mind, this pipeline was constructed with additional junction filtering tools inserted between pass 1 and pass 2 to minimize the incorporation of false positive novel juctions.

## Pipeline steps
1. [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. (optional) [fastp](https://github.com/OpenGene/fastp) Note: if fastp trimming is used, fastqc will run again on the trimmed fastq files
3. [STAR](https://github.com/alexdobin/STAR)
4. (2-pass mode only) [Portcullis](https://github.com/EI-CoreBioinformatics/portcullis) -- [documentation](https://portcullis.readthedocs.io/en/latest/)
5. (2-pass mode only) [annofilter_junctions.py](https://github.com/pfizer-rd/annofilter-junctions)
6. (2-pass mode only) STAR pass 2
7. [Salmon](https://github.com/COMBINE-lab/salmon) quant
8. (optional) [MAJIQ](https://majiq.biociphers.org/)
9. [Picard CollectRnaSeqMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard-)
10. [MultiQC](https://multiqc.info/)

## Usage
1. Create an input directory with a valid SampleSheet.csv (Illumina v1 or v2 acceptable) at the top level, along with a folder called `fastq` containing the fastq files to be processed.
2. Add entries to the header section of SampleSheet.csv to specify the strandedness (forward, reverse, unstranded), orientation (e.g., FR, RF, RR, FF, UNKNOWN), and read type (single, paired) of the reads. For example:
    ```
    [Header]

    Stranded,reverse
    Orientation,FR
    ReadType,paired
    ```
3. Provide required parameters. This is most easily accomplished using a params.yaml file:
    - `input`: [string] path to the input directory
    - `output`: [string] path to the desired output directory
    - `genome_id`: [string] the ID for the genome to use (see genome management section below)
    - `genome_registry`: [string] path to the genome registry to use for this run (see genome management section below)
    - `system_name`: [string] see genome management section below
4. (Optional) Provide optional pipeline parameters, again most easily done via a params.yaml file.
    - `allow_missing_samples`: Allow a run to continue with a warning if fastq files indicated in the SampleSheet.csv are not found (default: false)
    - `trim_reads`: [boolean] Use fastp to trim reads (default: false)
    - `star_2pass`: [boolean] Run STAR in 2-pass alignment mode (default: true)
    - `junctions`: [string] Path to an optional STAR junctions file to include in alignment
    - `read_length`: [integer] Used by STAR for genome indexing steps. As noted in STAR documentation, the default of 101 will work fine in the vast majority of cases. (default: 101)
    - `merge_bams`: [boolean] Merge all bam files from first pass alignment prior to analyzing newly discovered junctions. This can increase the data available for rare junctions, but can also produce a very large bam file depending on sequencing depth and the number of samples in the analysis. (default: true)
    - `rebuild_genome`: [boolean] Recreate the STAR genome index after the first pass alignment, incorporating newly discovered junctions. If false, the junctions will be added on-the-fly once for every sample. If true, this step will only be performed once and then the updated genome will be re-used for every sample. This can reduce computational resources used when there are more than a few samples being processed. (default: true)
    - `add_genes`: [string] ID of user-defined genes to add to the selected genome as an initial step before alignment. This is particularly useful when analyzing genetically modified samples or experiments involving delivery of exogenous gene products like gene therapy.
    - `run_majiq`: [boolean] Perform alternative splicing analysis with MAJIQ. (default: false)
    - `high_depth`: [boolean] Increases the computational resources requested for alignment steps. This is particularly useful for high sequencing depth libraries, e.g., > 75 million reads per sample. (default: false)
5. Submit a job request for the head node. Typically 2 cpus and 8 GB of memory will suffice. Request the appropriate queue based on expected wall-time of the run. The command to execute the job from within the job submission script is `nextflow run pfizer-rd/nf-rnaseq -params-file <params.yaml>`. 

## Default settings

=======
Star non-default parameters are as follows:

| Parameter | 2-pass, pass 1 | 2-pass, pass 2 | 1-pass | Comment |
| --------- | --------------- | -------------- | ------------- | --------- |
| alignIntronMin  | 20  | 20 | 20 | ENCODE |
| alignIntronMax | 1200000 | 1200000 | 1200000 |  This includes the maximum intron length in humans, 1160410. outSJfilterIntronMaxVsReadN will filter out unannotated junctions that create very long introns with few reads supporting them |
| alignSJoverhangMin | 8 | 999999 | 8 | 8 is ENCODE. 999999 is to prevent new junctions from being brought in in pass 2 |
| alignSJDBoverhangMin | 1 | 1 | 1 | ENCODE |
| alignMatesGapMax | 1200000 | 1200000 | 1200000 | This includes the maximum intron length in humans, 1160410 |
| outFilterMultimapNmax | 20 | 20 | 20 | ENCODE |
| outFilterMismatchNmax | 10 | 9999 | 9999 | 10 is default. 9999 is similar to ENCODE (999), but higher to make sure this variable doesn't have an affect and this is obvious to the user |
| outFilterMismatchNoverLmax | 0.3 | 0.04 | 0.04 | 0.3 is default for the first pass. 0.04 for the second pass and 1-pass is based on ENCODE outFilterMismatchNoverReadLmax=0.04, but we feel that it should be applied to the portion of the read that maps rather than the entire read length in case the actual mapped portion is much shorter than read length |

## Additional Usage Details:
### Genome repository
This pipeline makes use of [genome-manager](https://github.com/pfizer-opensource/genome-manager) and this is currently the only way to provide genome files. 

### Supplementation of genome with custom user genes
In some instances it is desirable to supplement the reference genome with custom user-defined genes. For instance, in the setting of gene therapy, one might treat cells with a specific gene therapy vector. Adding the vector sequence to the reference genome enables quantification of vector-derived reads. Such analyses are facilitated here by use of custom genome management software (insert github link). This dependency must be used separately to build a genome repository or, alternatively, one can simply provide the location of one that has already been built with the genome management tool.

### Containers
All pipeline steps require appropriate containers for reproducibility. Some are publicly available (e.g., quay.io) but internally developed containers are not provided and will need to be rebuilt. 

### License
This project is licensed under the MIT License; however, please ensure compliance with the licenses of any third-party software used within this pipeline.

## Authors
Rob Moccia, Melis Atalar Aksit, Abby Hill
