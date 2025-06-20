"""
Script Name: search_peptides_in_fasta.py
Description: Search peptides in protein FASTA.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import re,os,argparse
import numpy as np
from multiprocessing import Pool
from functools import partial

protein_seq_all_id_dict = {}
peptide_count_dict = {}
sample_list = []

def get_protein_id_to_seq_dict(input_file, header_delimiter=None, header_index=0):
    protein_id_to_seq_dict = {}
    for line in open(input_file, 'r'):
        if line[0]=='>':
            protein_id = line.rstrip('\n')[1:]
            if header_delimiter and header_index >= 0:
                protein_id = protein_id.split(header_delimiter)[header_index]
        else:
            protein_id_to_seq_dict[protein_id] = protein_id_to_seq_dict.get(protein_id, '') + line.rstrip('\n').upper()
    return protein_id_to_seq_dict


def get_peptide_count_dicts(input_file, exclude_list=None):
    peptide_count_dict = {}
    sample_list = []
    header_dict = {}
    abun_idx = 1
    for line in open(input_file, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
                if arr[i].startswith('total_'):
                    if abun_idx > 1:
                        continue
                    abun_idx = i
            sample_list = arr[abun_idx:]
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[header_dict['peptide']]
        peptide_count_dict[peptide_seq] = arr[abun_idx:]
    return peptide_count_dict,sample_list


def outf_table(peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_ends, output_file, append=False):
    if not append:
        outf = open(output_file, 'w')
        out_line = '#%s\t%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id', '\t'.join(sample_list))
        outf.write(out_line)
    else:
        outf = open(output_file, 'a')
    for peptide_seq in sorted(peptide_all_mapping_dict, key=lambda k:peptide_count_dict[k], reverse=True):
        tx_list = peptide_all_mapping_dict[peptide_seq]
        count_list = peptide_count_dict[peptide_seq]
        start_list = peptide_all_mapping_dict_starts[peptide_seq]
        end_list = peptide_all_mapping_dict_ends[peptide_seq]
        out_line = '%s\t%s\t%s\t%s\t%s\n'%(peptide_seq, ','.join(start_list), ','.join(end_list), ','.join(tx_list), '\t'.join(count_list))
        outf.write(out_line)
    outf.close()


def get_peptide_mappings(peptide_seq, protein_seq_all_id_dict):
    all_tx_list = []
    all_tx_start_list = []
    all_tx_end_list = []
    for protein_seq, all_seq_list in protein_seq_all_id_dict.items():
        for match in re.finditer(peptide_seq, protein_seq):
            start = match.start()+1 #Python mode
            end = match.end()
            for seq_id in all_seq_list: 
                all_tx_list.append(seq_id)
                all_tx_start_list.append(str(start))
                all_tx_end_list.append(str(end))
                #break
            #break
    if len(all_tx_list) > 0:
        return peptide_seq, 1, all_tx_list, all_tx_start_list, all_tx_end_list
    else:
        return peptide_seq, 0, None, None, None


def outf_table(results, output_file, peptide_count_dict, unmapped_peptide_set=set(), append=False):
    if output_file:
        if not append:
            outf = open(output_file, 'w')
            out_line = '#%s\t%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id', '\t'.join(sample_list))
            outf.write(out_line)
        else:
            outf = open(output_file, 'a')
    for result in results:
        peptide_seq, flag, all_tx_list, all_tx_start_list, all_tx_end_list = result
        if flag == 0:
            unmapped_peptide_set.add(peptide_seq)
        else:
            if output_file:
                count_list = peptide_count_dict[peptide_seq]
                out_line = '%s\t%s\t%s\t%s\t%s\n'%(peptide_seq, ','.join(all_tx_start_list), ','.join(all_tx_end_list), ','.join(all_tx_list), '\t'.join(count_list))
                outf.write(out_line)
    if output_file:
        outf.close()
    return unmapped_peptide_set


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search peptides in protein FASTA.')
    parser.add_argument('-pep', '--peptide_file', help='Input peptide table. Required.', required=True)
    parser.add_argument('-ref', '--ref_file', help='Input FASTA file of reference protein sequences. Required', required=True)
    parser.add_argument('-o2', '--output_unmapped', help='Output peptides not in reference FASTA', type=str)
    parser.add_argument('-o1', '--output_mapped', help='Output peptides in reference FASTA, mappings are added as additional columns.', type=str)
    parser.add_argument('--header_delimiter', help='Delimiter in the header if the header is too long to cut off, e.g., | in UniProt and GENCODE.', type=str)
    parser.add_argument('--header_index', help='If delimiter is given, the header is cut to list and item of given index is used as new header, default=0.', default=0, type=int)
    parser.add_argument('--threads', help='Number of threads to use, default=1.', default=1, type=int)
    args = parser.parse_args()

    ref_file = args.ref_file
    peptide_file = args.peptide_file
    if args.output_unmapped:
        output_file_unmapped = args.output_unmapped
    else:
        output_file_unmapped = None

    if args.output_mapped:
        output_file_mapped = args.output_mapped
    else:
        output_file_mapped = None

    if args.header_delimiter:
        header_delimiter = args.header_delimiter
    else:
        header_delimiter = None

    header_index = args.header_index

    threads = int(args.threads)
    pool = Pool(threads)

    protein_id_to_seq_dict = get_protein_id_to_seq_dict(ref_file, header_delimiter, header_index)
    protein_seq_all_id_dict = {}
    for seq_id, protein_seq in protein_id_to_seq_dict.items():
        protein_seq_all_id_dict.setdefault(protein_seq, []).append(seq_id)

    print('Entries in reference: %d'%(len(protein_id_to_seq_dict)))
    print('Unique proteins in reference: %d'%(len(protein_seq_all_id_dict)))

    peptide_count_dict,sample_list = get_peptide_count_dicts(peptide_file) #, exclude_list=set(novel_peptide_count_dict.keys()))
    print('Total peptides in input file: %d'%(len(peptide_count_dict)))

    processed_count = 0
    batch_size = 10000
    peptide_list = sorted(peptide_count_dict.keys(), key=lambda k:peptide_count_dict[k], reverse=True)
    unmapped_peptide_set = set()
    for i in range(0, len(peptide_list), batch_size):
        selected_peptide_list = peptide_list[i:i+batch_size]
        processed_count += len(selected_peptide_list)
        func=partial(get_peptide_mappings, protein_seq_all_id_dict=protein_seq_all_id_dict)
        results = pool.map(func, selected_peptide_list)
        if i == 0:
            append = False
        else:
            append = True
        unmapped_peptide_set = outf_table(results, output_file_mapped, peptide_count_dict, unmapped_peptide_set=unmapped_peptide_set, append=append)
        print('Processed %d peptides.'%(processed_count))

    if output_file_unmapped:
        outf = open(output_file_unmapped, 'w')
        for line in open(peptide_file, 'r'):
            if line.startswith('#'):
                outf.write(line)
                continue
            arr = line.rstrip('\n').split('\t')
            if arr[0] in unmapped_peptide_set:
                outf.write(line)
        outf.close()

    print('%s\t%s\t%s'%('#sample', 'mapped_peptides', 'unmapped_peptides'))
    print('%s\t%d\t%d'%(peptide_file, len(peptide_count_dict)-len(unmapped_peptide_set), len(unmapped_peptide_set)))
