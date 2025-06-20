"""
Script Name: select_peptides_by_fasta.py
Description: Filter peptide table by peptide sequences in FASTA file.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import re,os,argparse
import numpy as np

def get_seq_from_fasta(input_file):
    seq_dict = {}
    for line in open(input_file, 'r'):
        if line.startswith('>'):
            seq_id = line.rstrip('\n')[1:]
        else:
            seq_dict[seq_id] = seq_dict.get(seq_id, '') + line.rstrip('\n').upper()
    return seq_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter peptide table by peptide sequences in FASTA file.')
    parser.add_argument('-i', '--input_table', help='Input peptide table.', type=str, required=True)
    parser.add_argument('-f', '--input_fasta', help='Input FASTA of peptide sequences.', type=str)
    parser.add_argument('--retaining', help='Retaining peptide entries in the input table found in the given FASTA file. Set to False to discard peptides found in the FASTA. Default=True.', choices=['True', 'False'], default='True', type=str)
    parser.add_argument('-o','--output_table', help='Output peptide table', type=str)

    args = parser.parse_args()
    input_table = args.input_table
    input_fasta = args.input_fasta
    output_table = args.output_table
    if args.retaining == 'False':
        retaining = False
    else:
        retaining = True

    peptide_seq_dict = get_seq_from_fasta(input_fasta)
    peptide_seq_set = set(peptide_seq_dict.values())

    output = open(output_table, 'w')
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[0]
        if peptide_seq in peptide_seq_set:
            if retaining:
                output.write(line)
        else:
            if not retaining:
                output.write(line)

    output.close()

