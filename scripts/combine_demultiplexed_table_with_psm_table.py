import re,os,argparse
import numpy as np
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine demultiplexed intensities to PSM table')
    parser.add_argument('-d', '--demultiplexed_file', help='Input demultiplexed table with peptide intensities. Required.', required=True)
    parser.add_argument('-p', '--psm_file', help='Input PSM table with peptide mappings. Required.', required=True)
    parser.add_argument('-o', '--output_file', help='Output file with demultiplexed intensities and peptide mappings. Required', required=True)
    args = parser.parse_args()
    
    demultiplexed_file = args.demultiplexed_file
    psm_file = args.psm_file
    output_file = args.output_file

    header_dict = {}
    peptide_mappings_dict = {}
    for line in open(psm_file):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict) == 0:
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        peptide_seq = re.sub('\[UNIMOD:\d+\]', '', arr[4].split('.')[1]).upper()
        mapping_id_all = ';'.join([_.replace('Protein_', '') for _ in arr[5:]])
        peptide_mappings_dict[peptide_seq] = mapping_id_all

    output = open(output_file, 'w')
    for line in open(demultiplexed_file):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            output.write('#peptide\tprotein_index\t%s\n'%('\t'.join(arr[1:])))
            continue
        peptide_seq = arr[0]
        peptide_mapping = peptide_mappings_dict.get(peptide_seq, None)
        if not peptide_mapping:
            continue
        output.write('%s\t%s\t%s\n'%(peptide_seq, peptide_mapping, '\t'.join(arr[1:])))
    output.close()

