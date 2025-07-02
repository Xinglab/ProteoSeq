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
    parser.add_argument('-gf', '--gencode_fasta', help='Input FASTA of GENCODE translations. The header must be | separated for gene id and name annotation. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table with annotation retrieved from GTF. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_table = args.input_table
    collapsed_fasta = args.collapsed_fasta
    gencode_fasta = args.gencode_fasta
    output_table = args.output_table

    collapsed_protein_seq_dict = {}
    for line in open(collapsed_fasta, 'r'):
        if line.startswith('>'):
            protein_idx = line.rstrip('\n')[1:]
        else:
            collapsed_protein_seq_dict[protein_idx] = collapsed_protein_seq_dict.get(protein_idx, '') + line.rstrip('\n')

    gene_tx_protein_seq_dict = {}
    for line in open(gencode_fasta, 'r'):
        if line.startswith('>'):
            arr = line.rstrip('\n').split('|')
            gene_tx_id = (arr[2], arr[6], arr[1], arr[5])
        else:
            gene_tx_protein_seq_dict[gene_tx_id] = gene_tx_protein_seq_dict.get(gene_tx_id, '') + line.rstrip('\n')

    protein_seq_all_gene_tx_dict = {}
    for gene_tx_id, protein_seq in gene_tx_protein_seq_dict.items():
        protein_seq_all_gene_tx_dict.setdefault(protein_seq, set()).add(gene_tx_id)

    collapsed_protein_all_gene_tx_dict = {}
    for protein_idx, protein_seq in collapsed_protein_seq_dict.items():
        gene_tx_set = protein_seq_all_gene_tx_dict.get(protein_seq, set())
        if len(gene_tx_set) == 0:
            continue
        collapsed_protein_all_gene_tx_dict[protein_idx] = gene_tx_set

    annotated_protein_tx_dict = {}
    header_dict = {}
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = arr[header_dict['protein_idx']]
        transcript_id = arr[header_dict['transcript_id']]
        if transcript_id != '':
            annotated_protein_tx_dict.setdefault(protein_idx, set()).add(transcript_id)

    output = open(output_table, 'w')
    checked_protein_idx_set = set()
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = arr[header_dict['protein_idx']]
        if protein_idx in checked_protein_idx_set:
            output.write(line)
            continue
        checked_protein_idx_set.add(protein_idx)
        annotated_tx_set = annotated_protein_tx_dict.get(protein_idx, set())
        all_gene_tx_set = collapsed_protein_all_gene_tx_dict.get(protein_idx, set())
        for gene_tx in all_gene_tx_set:
            if gene_tx[2] in annotated_tx_set:
                continue
            output.write(f'{protein_idx}\t{gene_tx[2]}\tGENCODE\t{gene_tx[2]}\t{gene_tx[3]}\t{gene_tx[0]}\t{gene_tx[1]}\n')
        output.write(line)

    output.close()

