"""
Script Name: annotate_table_by_uniprot.py
Description: Annotate index table with information from UniProt file.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,sys,argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate index table with information from UniProt file.')
    parser.add_argument('-i', '--input_table', help='Input table protein index to original ids. Required', type=str, required=True)
    parser.add_argument('-f', '--collapsed_fasta', help='Input FASTA of final collapsed protein db. Required', type=str, required=True)
    parser.add_argument('-uf', '--uniprot_fasta', help='Input FASTA of UniProt. The header must still contain GN for gene name annotation. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table with annotation retrieved from GTF. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_table = args.input_table
    collapsed_fasta = args.collapsed_fasta
    uniprot_fasta = args.uniprot_fasta
    output_table = args.output_table

    collapsed_protein_seq_dict = {}
    for line in open(collapsed_fasta, 'r'):
        if line.startswith('>'):
            protein_idx = line.rstrip('\n')[1:]
        else:
            collapsed_protein_seq_dict[protein_idx] = collapsed_protein_seq_dict.get(protein_idx, '') + line.rstrip('\n')

    uniprot_gene_name_dict = {}
    uniprot_protein_seq_dict = {}
    for line in open(uniprot_fasta, 'r'):
        if line.startswith('>'):
            arr = line.rstrip('\n').split('|')
            gene_name = re.findall('GN=(\S*)', arr[2])
            if not gene_name:
                uniprot_id_gene = (arr[1], '')
                uniprot_gene_name_dict[arr[1]] = ''
            else:
                uniprot_id_gene = (arr[1],  gene_name[0])
                uniprot_gene_name_dict[arr[1]] = gene_name[0]
        else:
            uniprot_protein_seq_dict[uniprot_id_gene] = uniprot_protein_seq_dict.get(uniprot_id_gene, '') + line.rstrip('\n')

    protein_seq_all_gene_tx_dict = {}
    for uniprot_id_gene, protein_seq in uniprot_protein_seq_dict.items():
        protein_seq_all_gene_tx_dict.setdefault(protein_seq, set()).add(uniprot_id_gene)

    collapsed_protein_all_gene_tx_dict = {}
    for protein_idx, protein_seq in collapsed_protein_seq_dict.items():
        uniprot_id_gene_set = protein_seq_all_gene_tx_dict.get(protein_seq, set())
        if len(uniprot_id_gene_set) == 0:
            continue
        collapsed_protein_all_gene_tx_dict[protein_idx] = uniprot_id_gene_set

    annotated_uniprot_protein_dict = {}
    header_dict = {}
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = arr[header_dict['protein_idx']]
        seq_id = arr[header_dict['seq_id']]
        if arr[header_dict['ORF_type']] == 'UniProt' and seq_id != '':
            annotated_uniprot_protein_dict.setdefault(protein_idx, set()).add(seq_id)

    output = open(output_table, 'w')
    checked_seq_id_set = set()
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = arr[header_dict['protein_idx']]
        seq_id = arr[header_dict['seq_id']]
        if arr[header_dict['ORF_type']] == 'UniProt':
            arr[header_dict['gene_name']] = uniprot_gene_name_dict.get(seq_id, '')
        output.write('\t'.join(arr) + '\n')
        if protein_idx in checked_seq_id_set:
            continue
        checked_seq_id_set.add(protein_idx)
        annotated_uniprot_protein_set = annotated_uniprot_protein_dict.get(protein_idx, set())
        all_uniprot_id_gene_set = collapsed_protein_all_gene_tx_dict.get(protein_idx, set())
        for uniprot_id_gene in all_uniprot_id_gene_set:
            if uniprot_id_gene[0] in annotated_uniprot_protein_set:
                continue
            output.write(f'{protein_idx}\t{uniprot_id_gene[0]}\tUniProt\t\t\t\t{uniprot_id_gene[1]}\n')

    output.close()

