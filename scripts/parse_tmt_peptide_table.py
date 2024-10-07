import re,os,argparse
import numpy as np
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse demultiplexed files.')
    parser.add_argument('-i', '--input_file', help='Input tsv table output from TextExporter of consensusXML file. Required.', required=True)
    parser.add_argument('-o', '--output_file', help='Output file with demultiplexed intensities of peptides. Required', required=True)
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    tmt_to_sample_dict = {}

    map_to_tmt_id_dict = {}
    for line in open(input_file):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[0] == 'MAP':
            map_to_tmt_id_dict[arr[1]] = arr[3]

    feature_n_to_peptide_dict = {}
    feature_n = 0
    tmt_int_col_dict = {}
    for line in open(input_file):
        arr = line.rstrip('\n').split('\t')
        if arr[0] == '#CONSENSUS':
            for i in range(len(arr)):
                matches = re.findall('intensity_(\d+)', arr[i])
                if len(matches) > 0:
                    tmt_id = map_to_tmt_id_dict[matches[0]]
                    tmt_int_col_dict[tmt_id] = i
            continue
        if arr[0] == 'CONSENSUS':
            feature_n += 1
            continue
        if arr[0] == 'PEPTIDE':
            peptide_seq = arr[5].strip('.')
            peptide_seq = peptide_seq.replace('(TMT6plex)', '')
            peptide_seq = peptide_seq.replace('(Carbamidomethyl)', '')
            peptide_seq = peptide_seq.replace('(Oxidation)', '')
            feature_n_to_peptide_dict[feature_n] = peptide_seq

    feature_n = 0
    peptide_tmt_int_dict = {}
    for line in open(input_file):
        arr = line.rstrip('\n').split('\t')
        if arr[0] == 'CONSENSUS':
            feature_n += 1
            peptide_seq = feature_n_to_peptide_dict.get(feature_n, None)
            if not peptide_seq:
                continue
            for tmt_id, col_id in tmt_int_col_dict.items():
                intensity = float(arr[col_id])
                peptide_tmt_int_dict[peptide_seq][tmt_id] = peptide_tmt_int_dict.setdefault(peptide_seq, {}).get(tmt_id, 0) + intensity

    output = open(output_file, 'w')
    tmt_list = sorted(tmt_int_col_dict.keys())
    if tmt_to_sample_dict:
        sample_list = [tmt_to_sample_dict.get(_, _) for _ in tmt_list]
    else:
        sample_list = tmt_list[:]

    output.write('#peptide\ttotal_intensity\t%s\n'%('\t'.join(sample_list)))
    for peptide_seq in sorted(peptide_tmt_int_dict, key=lambda k:sum(peptide_tmt_int_dict[k].values()), reverse=True):
        tmt_int_dict = peptide_tmt_int_dict[peptide_seq]
        int_list = [tmt_int_dict.get(_, 0) for _ in tmt_list]
        total_int = sum(int_list)
        int_list_str = ['%.2f'%(_) for _ in int_list]
        output.write('%s\t%.2f\t%s\n'%(peptide_seq, total_int, '\t'.join(int_list_str)))

    output.close()
