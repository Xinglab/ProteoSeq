import os,argparse,re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collapse multiple PSM tables.')
    parser.add_argument('-i', '--input_file_list', help='List of input PSM tables, space splitted. Required', nargs='+', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=str)
    parser.add_argument('--sample_list', help='Input the sample name lists for the input files. PSMs in each input will be collapsed. If not given, all columns in the input table are output', nargs='+', type=str)
    parser.add_argument('--combine_sample', help='Combine the runs in each sample file. Default=True', default='True', choices=['True', 'False'], type=str)
    args = parser.parse_args()
    
    input_file_list = args.input_file_list
    output_file = args.output_file
    if args.combine_sample == 'True':
        combine_sample = True
    else:
        combine_sample = False
    if args.sample_list:
        if len(args.sample_list) == len(input_file_list):
            sample_list = args.sample_list
        else:
            exit('Length of input file list is not equal to the length of sample list.')
    else:
        sample_list = [_.split('/')[-1].replace('.sig_psms_count.txt', '') for _ in input_file_list]

    peptide_all_mapping_dict = {}
    all_samples_psm_count_dict = {}
    all_ms_run_list = []
    for input_i in range(len(input_file_list)):
        input_file = input_file_list[input_i]
        sample = sample_list[input_i]
        psm_count_dict = {}
        header_dict = {}
        ms_run_list = []
        count_column = -1
        count_int = True
        for line in open(input_file, 'r'):
            arr = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                for i in range(len(arr)):
                    header_dict[arr[i]] = i
                    if arr[i] not in ['#peptide', 'peptide', 'protein_index', 'total_count', 'total_intensity', 'starts', 'ends', 'tx_id']:
                        ms_run_list.append(arr[i])
                        all_ms_run_list.append(arr[i])
                    elif arr[i] == 'total_count':
                         count_column = i
                         count_int = True
                    elif arr[i] == 'total_intensity':
                         count_column = i
                         count_int = False
                continue
            peptide_seq = arr[0]
            peptide_all_mapping_dict.setdefault(peptide_seq, set()).add(arr[1])
            if combine_sample:
                if count_int:
                    all_samples_psm_count_dict[peptide_seq][sample] = all_samples_psm_count_dict.setdefault(peptide_seq, {}).get(sample, 0) + int(arr[count_column])
                else:
                    all_samples_psm_count_dict[peptide_seq][sample] = all_samples_psm_count_dict.setdefault(peptide_seq, {}).get(sample, 0) + float(arr[count_column])
            else:
                for ms_run in ms_run_list:
                    if count_int:
                        all_samples_psm_count_dict.setdefault(peptide_seq, {})[ms_run] = int(arr[header_dict[ms_run]])
                    else:
                        all_samples_psm_count_dict.setdefault(peptide_seq, {})[ms_run] = float(arr[header_dict[ms_run]])

    peptide_all_mapping_combined_dict = {}
    for peptide_seq,mapping_set in peptide_all_mapping_dict.items():
        if len(mapping_set)>1:
            print('%s has multiple annotations:%s'%(peptide_seq, ' '.join(mapping_set)))
        for mapping in mapping_set:
            for m in mapping.split(';'):
                protein_idx = int(m.replace('Protein_', ''))
                peptide_all_mapping_combined_dict.setdefault(peptide_seq, set()).add(protein_idx)

    output = open(output_file, 'w')
    if combine_sample:
        output_sample_list = sample_list[:]
    else:
        output_sample_list = all_ms_run_list[:]

    if count_int:
        output.write('#peptide\tprotein_index\ttotal_count\t%s\n'%('\t'.join(output_sample_list)))
    else:
        output.write('#peptide\tprotein_index\ttotal_intensity\t%s\n'%('\t'.join(output_sample_list)))
    for peptide_seq in sorted(all_samples_psm_count_dict, key=lambda k:sum(all_samples_psm_count_dict[k].values()), reverse=True):
        mapping = ';'.join([str(_) if type(_)==int else _ for _ in sorted(peptide_all_mapping_combined_dict[peptide_seq])])
        count_list = [all_samples_psm_count_dict[peptide_seq].get(sample, 0) for sample in output_sample_list]
        total_count = sum(count_list)
        if count_int:
            count_list = ['%d'%_ for _ in count_list]
        else:
            count_list = ['%.2f'%_ for _ in count_list]
        if count_int:
            output.write('%s\t%s\t%d\t%s\n'%(peptide_seq, mapping, total_count, '\t'.join(count_list)))
        else:
            output.write('%s\t%s\t%.2f\t%s\n'%(peptide_seq, mapping, total_count, '\t'.join(count_list)))

    output.close()


