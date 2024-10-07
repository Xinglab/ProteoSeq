import os,argparse,re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collapse PIN outputs to tabular table.')
    parser.add_argument('-i', '--input_file_list', help='List of input pin files with PSMs, space splitted. Required', nargs='+', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=str)
    args = parser.parse_args()
    
    input_file_list = args.input_file_list
    output_file = args.output_file

    peptide_tx_dict = {}
    peptide_sample_count_dict = {}
    sample_set = set()
    for input_file in input_file_list:
        sample = input_file.split('/')[-1].replace('.sig.pin', '')
        sample = sample.replace('.pin', '')
        sample_set.add(sample)
        header_dict = {}
        for line in open(input_file, 'r'):
            arr = line.rstrip('\n').split('\t')
            if len(header_dict) == 0:
                for i in range(len(arr)):
                    header_dict[arr[i]] = i
                continue
            peptide_seq = arr[header_dict['Stripped.Sequence']]
            peptide_sample_count_dict[peptide_seq][sample] = peptide_sample_count_dict.setdefault(peptide_seq, {}).get(sample, 0) + 1
            tx_id_all = arr[header_dict['Protein.Group']]
            peptide_tx_dict.setdefault(peptide_seq, set()).add(tx_id_all)

    peptide_to_protein_dict = {}
    for peptide_seq,mapping_set in peptide_tx_dict.items():
        if len(mapping_set)>1:
            print('%s has multiple annotations:%s'%(peptide_seq, ' '.join(mapping_set)))
        for mapping in mapping_set:
            for m in mapping.split(';'):
                protein_idx = int(m.replace('Protein_', ''))
                peptide_to_protein_dict.setdefault(peptide_seq, set()).add(protein_idx)

    output = open(output_file, 'w')
    sample_list = sorted(sample_set)
    output.write('#peptide\tprotein_index\ttotal_count\t%s\n'%('\t'.join(sample_list)))
    for peptide_seq in sorted(peptide_sample_count_dict, key=lambda k:sum(peptide_sample_count_dict[k].values()), reverse=True):
        count_list = [peptide_sample_count_dict[peptide_seq].get(sample, 0) for sample in sample_list]
        total_count = sum(count_list)
        count_list = [str(_) for _ in count_list]
        protein_idx_output = [str(_) if type(_)==int else _ for _ in sorted(peptide_to_protein_dict[peptide_seq])]
        output.write('%s\t%s\t%d\t%s\n'%(peptide_seq, ';'.join(protein_idx_output), total_count, '\t'.join(count_list)))        

    output.close()


