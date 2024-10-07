import re,os,argparse
import numpy as np


def get_peptide_dict(input_file):
    sample_len_count_dict = {}
    sample_pep_set_dict = {}
    header_dict = {}
    header_rev_dict = {}
    sample_list = []
    abun_idx = 1
    for line in open(input_file, 'r'):
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
                header_rev_dict[i] = arr[i]
                if arr[i].startswith('total_'):
                    if abun_idx > 1:
                        continue
                    abun_idx = i
            sample_list = arr[abun_idx+1:]
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[0]
        peptide_len = len(peptide_seq)
        for i in range(abun_idx+1, len(arr)):
            peptide_abun = float(arr[i])
            sample = header_rev_dict[i]
            if peptide_abun > 0:
                sample_pep_set_dict.setdefault(sample, set()).add(peptide_seq)
                sample_len_count_dict[sample][peptide_len] = sample_len_count_dict.setdefault(sample, {}).get(peptide_len, 0) + peptide_abun
    return sample_len_count_dict,sample_list,sample_pep_set_dict

def merge_dicts(d1, d2):
    new_d = d1.copy()
    for k, v in d2.items():
        if isinstance(v, dict):
            for k2,v2 in v.items():
                new_d[k][k2] = new_d.setdefault(k, {}).get(k2, 0) + v2
        elif isinstance(v, set):
            new_d[k] = new_d.setdefault(k, set()) | v
        else:
            print('Wrong input dictionaries.')
            return False
    return new_d

def get_len_abun_dicts(sample_len_count_dict):
    total_len_dict = {}
    total_abun_dict = {}
    for sample, peptide_abun_dict in sample_len_count_dict.items():
        for peptide_len, peptide_abun in peptide_abun_dict.items():
            total_len_dict[sample] = total_len_dict.get(sample, 0) + peptide_len*peptide_abun
            total_abun_dict[sample] = total_abun_dict.get(sample, 0) + peptide_abun
    return total_len_dict,total_abun_dict

def get_pep_out_line(sample_pep_set_dict, sample_list):
    pep_count_list = []
    all_pep_set = set()
    for sample in sample_list:
        pep_set = sample_pep_set_dict.get(sample, set())
        pep_count_list.append(str(len(pep_set)))
        all_pep_set |= pep_set
    return '%d\t%s'%(len(all_pep_set), '\t'.join(pep_count_list))


def get_peptide_abun_out_line(total_abun_dict, sample_list):
    abun_list = []
    total_abun = 0
    for sample in sample_list:
        abun = total_abun_dict.get(sample, 0)
        total_abun += abun
        abun_list.append('%d'%(abun))
    return '%d\t%s'%(total_abun, '\t'.join(abun_list))


def get_length_out_line(total_len_dict, total_abun_dict, sample_list):
    if sum(total_abun_dict.values()) == 0:
        return '0\t%s'%('\t'.join(['0']*len(sample_list)))
    total_len, total_abun = 0, 0
    avg_len_list = []
    for sample in sample_list:
        sample_len = total_len_dict.get(sample, 0)
        abun = total_abun_dict.get(sample, 0)
        total_len += sample_len
        total_abun += abun
        if abun == 0:
            avg_len_list.append('0')
        else:
            avg_len_list.append('%.2f'%(sample_len/abun))
    if total_abun == 0:
        return '0\t%s'%('\t'.join(avg_len_list))
    return '%.2f\t%s'%(total_len/total_abun, '\t'.join(avg_len_list))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize peptides.')
    parser.add_argument('-i', '--input_tables', help='Input peptide tables, space splitted.', nargs='+', type=str, required=True)
    parser.add_argument('-n', '--input_names', help='Input names for the tables, space splitted.', nargs='+', type=str, required=True)
    parser.add_argument('-o','--output_file', help='Output summary reports', type=str)

    args = parser.parse_args()
    
    input_tables_list = args.input_tables
    input_names_list = args.input_names
    output_file = args.output_file

    if len(input_tables_list) != len(input_names_list):
        exit('Input tables and names do not have same length.')

    sample_len_count_dict_all = {}
    sample_list_all = {}
    sample_pep_set_dict_all = {}
    for i in range(len(input_tables_list)):
        input_table = input_tables_list[i]
        input_name = input_names_list[i]
        sample_len_count_dict, sample_list, sample_pep_set_dict = get_peptide_dict(input_table)
        if not sample_len_count_dict_all.get(input_name, None):
            sample_len_count_dict_all[input_name] = sample_len_count_dict.copy()
        else:
            sample_len_count_dict_all[input_name] = merge_dicts(sample_len_count_dict_all[input_name], sample_len_count_dict)
        if not sample_pep_set_dict_all.get(input_name, None):
            sample_pep_set_dict_all[input_name] = sample_pep_set_dict.copy()
        else:
            sample_pep_set_dict_all[input_name] = merge_dicts(sample_pep_set_dict_all[input_name], sample_pep_set_dict)
        sample_list_all[input_name] = sorted(set(sample_list_all.get(input_name, []))|set(sample_list))

    sample_set = set()
    for input_name, input_sample_list in sample_list_all.items():
        sample_set |= set(input_sample_list)
    del sample_list_all
    sample_list = sorted(sample_set)
    del sample_set

    output = open(output_file, 'w')
    output.write('## Number of samples: %d\n'%(len(sample_list)))
    output.write('\n')
    output.write('## Number of peptides in all/each sample:\n')
    output.write('#Peptides_types\tPeptides_count\t%s\n'%('\t'.join(sample_list)))
    for input_name in set(input_names_list):
        out_line = get_pep_out_line(sample_pep_set_dict_all[input_name], sample_list)
        output.write('%s\t%s\n'%(input_name, out_line))

    del sample_pep_set_dict_all
    print('Finished writing peptide summary to the report.')

    total_len_dict_all = {}
    total_abun_dict_all = {}
    for input_name in input_names_list:
        total_len_dict_all[input_name], total_abun_dict_all[input_name] = get_len_abun_dicts(sample_len_count_dict_all[input_name])

    del sample_len_count_dict_all

    output.write('\n')
    output.write('## Number of PSMs in all/each sample:\n')
    output.write('#Peptides_types\tTotal_PSM\t%s\n'%('\t'.join(sample_list)))
    for input_name in set(input_names_list):
        out_line = get_peptide_abun_out_line(total_abun_dict_all[input_name], sample_list)
        output.write('%s\t%s\n'%(input_name, out_line))

    print('Finished writing PSM summary to the report.')

    output.write('\n')
    output.write('## Average peptide length in all/each sample:\n')
    output.write('#Peptides_types\tAvg_len\t%s\n'%('\t'.join(sample_list)))
    for input_name in set(input_names_list):
        out_line = get_length_out_line(total_len_dict_all[input_name], total_abun_dict_all[input_name], sample_list)
        output.write('%s\t%s\n'%(input_name, out_line))

    output.close()

