import csv
# import re
import argparse

def parse_majiq( deltapsi_file, output_file ):
    # Open the TSV file
    with open(deltapsi_file, 'r') as tsvfile:
        # Create a TSV reader object
        reader = csv.reader(tsvfile, delimiter='\t')
   
        # Create a new TSV file to write the results
        with open(output_file, 'w', newline='') as outfile:
            # Create a TSV writer object
            writer = csv.writer(outfile, delimiter='\t',quoting=csv.QUOTE_ALL)

            # copy and add to header
            header_row = next(reader)
            header_row.extend(list(['LSV_has_ir','ref_exon_type','ref_exon_SJ_num','alt_exon_info','alt_exon_num','exons_skipped']))
            writer.writerow(header_row)

            # Iterate over each remaining row in the TSV file
            for row in reader:
                # print(row)
                # Get the number of semicolons in the first column of the row
                num_events = row[3].count(';') +1
                # print('current LSV number of events:',num_events)

                # go through the value in each column of that row
                for curr_event in range(0,num_events):
                    # make a new row identical to the original one
                    curr_event_row = row.copy()
                    for curr_col, values in enumerate(row):
                        # if column is not empty on that row
                        if values:
                            curr_lsv_has_ir = True
                            # Count how many values there should be in each column
                            num_event_values = values.count(';') +1
                            if (num_event_values >1):
                                curr_junction_values = values.split(';')
                                curr_event_row[curr_col] = curr_junction_values[curr_event]
                                # print(curr_junction_values[curr_event])
                            elif (values.count('|') > 1):
                                # this happens if there are multiple events
                                curr_lsv_type_info = values.split('|')
                                curr_lsv_event_type = curr_lsv_type_info[curr_event+1]
                                curr_event_row[curr_col] = curr_lsv_event_type
                                # print(curr_lsv_event_type)
                                if (curr_lsv_event_type == 'i'):
                                    curr_alt_exon_info = 'i'
                                    curr_ref_exon_SJ_num = '0'
                                    curr_alt_exon_num = '0'
                                    curr_exons_skipped = '-1'
                                else:
                                    # print(curr_lsv_event_type)
                                    curr_alt_exon_info = curr_lsv_event_type.split('e')[1]
                                    curr_ref_exon_SJ_num = curr_lsv_event_type.split('e')[0]
                                    curr_alt_exon_num = curr_alt_exon_info.split('.')[0]
                                    if (curr_lsv_type_info[0] == 's'):
                                        # print(curr_alt_exon_num)
                                        curr_exons_skipped = str(int(curr_alt_exon_num)-1)
                                    elif (curr_lsv_type_info[0] == 't'):
                                        curr_exons_skipped = str(int(row[9])-1-int(curr_alt_exon_num))
                                        # print(row[8],row[9],curr_exons_skipped)
                                    else:
                                        curr_exons_skipped = 'na'

                            elif (values.count('|') == 1):
                                # this happens if there are too many events (>30?) and it just says 'na'
                                # make all the additional columns 'na'
                                curr_lsv_type_info = values.split('|')
                                curr_lsv_event_type = curr_lsv_type_info[1]
                                curr_event_row[curr_col] = curr_lsv_type_info[1]
                                curr_alt_exon_info = curr_lsv_type_info[1]
                                curr_ref_exon_SJ_num = curr_lsv_type_info[1]
                                curr_alt_exon_num = curr_lsv_type_info[1]
                                curr_exons_skipped = curr_lsv_type_info[1]
                            else:
                                curr_event_row[curr_col] = row[curr_col]
                        else:
                            # empty columns mean intron retention coordinates are absent
                            curr_lsv_has_ir = False
                    # reminder of the order: list(['LSV_has_ir','ref_exon_type','ref_exon_SJ_num','alt_exon_info','alt_exon_num','exons_skipped']))
                    curr_event_row.extend(list([curr_lsv_has_ir,curr_lsv_type_info[0]]))
                    curr_event_row.extend(list([curr_ref_exon_SJ_num,curr_alt_exon_info,curr_alt_exon_num,curr_exons_skipped]))

                    # Write the new row to the new TSV file
                    # print(curr_event_row)
                    writer.writerow(curr_event_row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser( 
        description = 'Script to convert MAJIQ deltapsi.tsv output from 1 row per LSV to 1 row per splicing event (splice junction or intron retention event).')
    parser.add_argument( '--deltapsi', '-d', help = 'Full path to the contrast.deltapsi.tsv file (NOT .voila.tsv)' )
    parser.add_argument( '--output', '-o', help = 'Full path to write output file, usually contrast.deltapsi.split.tsv' )
    
    args = parser.parse_args()
    
    parse_majiq( args.deltapsi, args.output )