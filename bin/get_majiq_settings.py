#!/usr/bin/env python3

__author__ = 'Abby Hill'
__version__ = '0.1'

import argparse
import os.path
from collections import defaultdict
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Script to convert a SampleSheetNextflow.csv from parse_sample_sheet.py containing a column majiq_group into a majiq_settings.ini for running MAJIQ build using Nextflow pipeline.')
    parser.add_argument('--samplesheetnextflow',
        help = 'reformatted SampleSheet.csv from PARSE_SAMPLESHEET')
    parser.add_argument('--majiqcontrasts',
        help = 'majiq_contrasts.csv provided by user')
    parser.add_argument('--genome',
        help = 'name of genome to use for links to UCSC genome browser, may differ from genome manager genome names')
    parser.add_argument('--outdir', '-o',
        help = 'directory to write modified majiq_settings.ini')
    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    if args.outdir:
        outfile_settings = os.path.join(args.outdir, 'majiq_settings.ini')
        outfile_voilacommands = os.path.join(args.outdir, 'voila_view_cmds.txt')
    else:
        outfile_settings = 'majiq_settings.ini'
        outfile_voilacommands = 'voila_view_cmds.txt'

    majiq_groups = defaultdict(list)
    stranded_values = []
    with open(args.samplesheetnextflow, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            majiq_group = row['majiq_group']
            sample_name = row['sample_name']
            # Append the sample ID to the list of sample IDs for that majiq_group
            majiq_groups[majiq_group].append(sample_name+'.Aligned.sortedByCoord.out')
            stranded_values.append(row['stranded'])

    if stranded_values.count(stranded_values[0]) == len(stranded_values):
        stranded = stranded_values[0]
    else:
        raise ValueError("Error: Not all values in 'stranded' are identical.")
    if stranded == 'unstranded':
        stranded = 'none'
    valid_stranded_values = set(['forward', 'reverse', 'none'])
    if stranded not in valid_stranded_values:
        raise ValueError("Error: Invalid value for stranded. MAJIQ strandness must be 'forward', 'reverse', or 'none'.")


    # Open the settings.ini file for writing
    with open(outfile_settings, 'w') as ini_file:
        ini_file.write("[info]\ngenome={}\nbamdirs=./\n".format(args.genome))
        ini_file.write("strandness={}\n".format(stranded))
        ini_file.write("[experiments]\n")
        # Iterate over the majiq_groups and their sample IDs
        for majiq_group, sample_names in majiq_groups.items():
            # Write the majiq_group and its sample IDs to the settings.ini file
            ini_file.write("{}={}\n".format(majiq_group,','.join(sample_names)))

    majiq_contrasts = defaultdict(list)
    with open(args.majiqcontrasts, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            contrast_name = row['contrast_name']
            majiq_contrasts[contrast_name].append(row['reference'])
            majiq_contrasts[contrast_name].append(row['experiment'])


   # Open the voila.bsub file for writing
    with open(outfile_voilacommands, 'w') as voila_file:
        voila_file.write("module load ib python/3.8.0\n")
        voila_file.write("source /hpc/grid/shared/ngsapp/majiq-commercial/bin/activate\n")
        voila_file.write("port_number=2021\n")
        # Iterate over the majiq_contrasts
        i = 0
        for names, conditions in majiq_contrasts.items():
            # Write the majiq_group and its sample IDs to the settings.ini file
            voila_file.write("bsub -q long voila view build_out/splicegraph.sql  deltapsi_out/{}.deltapsi.voila --num-web-workers 2 --host 0.0.0.0 --port  $((port_number + {}))\n".format('-'.join(conditions),i))
            i+=1

    # test_samplesheetnextflow = '/lustre/workspace/home/hilla34/Projects/RBM20/Fang_Lu_RBM20_23mRNAseq/pre_run/processed_nextflow_majiq/metadata/SampleSheetNextflow.csv'
    # test_majiqcontrasts = '/lustre/workspace/home/hilla34/Projects/RBM20/Fang_Lu_RBM20_23mRNAseq/pre_run/majiq_contrasts.csv'