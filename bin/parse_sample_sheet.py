#!/usr/bin/env python3

__author__ = 'Rob Moccia'
__version__ = '0.1'

import argparse
import sys
from itertools import groupby
import os.path

class ValueMissingError(Exception):
    pass

class ReadTypeOrientationMismatchError(Exception):
    pass


def generate_minimal_samplesheet(sample_ids, single_end, stranded):
    '''
    Create a minimal sample sheet sufficient for Nextflow run if no sample sheet provided.
    '''
    colnames = ['sample_id', 'single_end', 'stranded']
    yield ','.join(colnames)
    for sample in sample_ids:
        yield f'{sample},{str(single_end).lower()},{stranded}'

def parse_illumina_samplesheet(samplesheet, stranded = None, read_type = None, orientation = 'UNKNOWN'):
    '''
    Parse a SampleSheet and modify it to contain columns required for Nextflow run.
    Expects Illumina format and optional comma separated key-value pairs in [Header]:
        Stranded, <'unstranded', 'forward', 'reverse'>
        ReadType, <'SE', 'single', 'single-end', 'PE', 'paired', 'paired-end'>
        Orientation, <>
    '''
    with open(samplesheet, 'r') as f:
        # section headers are enclosed by square brackets so group the file
        # into sections based on this and store them in a dictionary
        groups = groupby(f,
            key = lambda x: x.strip().startswith('['))
        sections = dict()
        for k, v in groups:
            if k:
                section_name = next(v).strip('[]\n, ')
                section_data = [row.strip() for row in next(groups)[1] if ',' in row]
                sections[section_name] = section_data

        # header should contain comma delimited key-value pairs
        # parse these and store in a dictionary
        header_dict = dict()
        for line in sections['Header']:
            key, value = line.split(',')[:2]
            header_dict[key] = value.strip('\'"')

        # for key fields like Stranded, the SampleSheet takes precedence
        # if it's not declared in SampleSheet header, check if it was supplied
        # as a function argument
        stranded_field = header_dict.get('Stranded')
        valid_stranded = ['unstranded', 'forward', 'reverse', 'none']
        if stranded_field is None:
            if stranded is None:
                raise ValueMissingError("If 'Stranded' is not specified in SampleSheet header it must be supplied as an argument.")
            else:
                stranded_field = stranded.lower()
        else:            
            print("Using 'Stranded' specified in SampleSheet")
            if stranded is not None:
                print(f'Script called with argument "stranded" = {stranded}. Ignoring argument and using value provided in SampleSheet.')
            stranded_field = stranded_field.lower()
        if stranded_field not in valid_stranded:
            raise ValueError(f'Stranded must be one of {valid_stranded}')
        # allow for the user to say "none" but convert to the pipeline convention of "unstranded"
        if stranded_field == 'none':
            stranded_field = 'unstranded'

        # ReadType
        read_type_field = header_dict.get('ReadType')
        valid_read_types = {'se': 'single', 'single': 'single', 'single-end': 'single',
            'pe': 'paired', 'paired': 'paired', 'paired-end': 'paired'}
        if read_type_field is None:
            if read_type is None:
                raise ValueMissingError("If 'ReadType' is not specified in SampleSheet header it must be supplied as an argument.")
            else:
                read_type_field = read_type.lower()
        else:
            print("Using 'ReadType' specified in SampleSheet.")
            if read_type is not None:
                print(f'Script called with argument "read_type" = {read_type}. Ignoring argument and using value provided in SampleSheet.')
            read_type_field = read_type_field.lower()
        if read_type_field not in valid_read_types.keys():
            raise ValueError(f'ReadType must be one of {valid_read_types.keys()}')
        # standardize vocabulary
        read_type_field = valid_read_types[read_type_field]
        # convert to boolean
        is_paired_end = read_type_field == 'paired'

        # Orientation is used by Portcullis but also be 'UNKNOWN'
        valid_se_orientations = ['F', 'R', 'UNKNOWN']
        valid_pe_orientations = ['FR', 'RF', 'FF', 'RR', 'UNKNOWN']
        valid_orientations = set(valid_se_orientations).union(set(valid_pe_orientations))
        orientation_field = header_dict.get('Orientation')
        if orientation_field is None:
            print(f"'Orientation' not specified in SampleSheet header. Using argument 'orientation' = {orientation}")
            orientation_field = orientation
        else:
            print("Using 'Orientation' specified in SampleSheet.")
            orientation_field = orientation_field.upper()
        if orientation_field not in valid_orientations:
            raise ValueError(f'Orientation must be one of {valid_orientations}')
        # check for incompatible ReadType and Orientation
        if is_paired_end and not orientation_field in valid_pe_orientations \
              or (not is_paired_end and not orientation_field in valid_se_orientations):
            raise ReadTypeOrientationMismatchError(f'ReadType {read_type_field} and Orientation {orientation_field} are incompatible with each other.')

        # check for Illumina version 2 file format in which case data will be in [BCLConvert_Data] rather than [Data]
        if header_dict.get('FileFormatVersion') == '2':
            print('Detected Illumina SampleSheet file-format version 2.')
            read_data = sections['BCLConvert_Data']
        else:
            try:
                read_data = sections['Data']
            except:
                print(header_dict)
                raise

        # first entry in Data section will be the column names
        colnames = read_data[0].split(',') + ['read_type', 'stranded', 'orientation']
        colnames = [col.lower() for col in colnames]
        # Sample_ID is not required to be the first column in the SampleSheet
        # but it will be required to be the first column in the output
        sample_col = colnames.index('sample_id')
        output_order = list(range(0, len(colnames)))
        if sample_col != 0:
            output_order.insert(0, output_order.pop(sample_col))
            colnames = [colnames[i] for i in output_order]

        yield ','.join(colnames)

        # parse the data rows
        for sample in read_data[1:]:
            if sample != '':
                fields = sample.split(',') + [read_type_field, stranded_field, orientation_field]
                fields = [fields[i].strip('\'"') for i in output_order]
                yield ','.join(fields)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Script to convert and Illumina SampleSheet.csv into a format suitable for Nextflow pipeline.')
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('--samplesheet',
        help = 'Illumina format SampleSheet.csv')
    group.add_argument('--samples',
        help = 'list of sample names')
    parser.add_argument('--outdir', '-o',
        help = 'directory to write modified SampleSheet.csv')
    parser.add_argument('--stranded', '-s',
        choices = ['unstranded', 'forward', 'reverse', None],
        required = '--samples' in sys.argv,
        default = None,
        help = "read stranded")
    parser.add_argument('--read_type',
        choices = ['single', 'paired'],
        required = '--samples' in sys.argv,
        default = None,
        help = 'single or paired reads')
    parser.add_argument('--orientation',
        choices = ['F', 'R', 'FR', 'RF', 'FF', 'RR', 'UNKNOWN'],
        required = '--samples' in sys.argv,
        default = 'UNKNOWN',
        help = 'orientation of reads')
    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()

    if args.outdir:
        outfile = os.path.join(args.outdir, 'SampleSheetNextflow.csv')
    else:
        outfile = 'SampleSheetNextflow.csv'

    with open(outfile, 'w') as f:
        for line in parse_illumina_samplesheet(
            args.samplesheet,
            stranded = args.stranded,
            read_type = args.read_type):
            f.write(line + '\n')

    # test_samplesheet = '/home/moccir/nextflow_test/SampleSheetFullSRA.csv'
    # sample_ids = ['SRR101', 'SRR102', 'SRR103']
    # for line in generate_minimal_samplesheet(sample_ids, False, 'none'):
    #     print(line)
