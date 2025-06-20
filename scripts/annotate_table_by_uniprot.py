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
    parser.add_argument('-f', '--uniprot_fasta', help='Input FASTA of UniProt. The header must still contain GN for gene name annotation. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table with annotation retrieved from GTF. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_table = args.input_table
    uniprot_fasta = args.uniprot_fasta
    output_table = args.output_table

    gene_name_dict = {}
    for line in open(uniprot_fasta, 'r'):
        if line.startswith('>'):
            arr = line.rstrip('\n').split('|')
            seq_id = arr[1]
            gene_name = re.findall('GN=(\S*)', arr[2])
            if len(gene_name) == 0:
                continue
            gene_name_dict[seq_id] = gene_name[0]

    output = open(output_table, 'w')
    header_dict = {}
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            if header_dict.get('gene_name', -1) < 0:
                output.write('%s\t%s\n'%(line.rstrip('\n'), 'gene_name'))
            else:
                output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        seq_id = arr[header_dict['seq_id']]
        gene_name = gene_name_dict.get(seq_id, '')
        gene_name_column = header_dict.get('gene_name', -1)
        if gene_name_column < 0:
            output.write('%s\t%s\n'%(line.rstrip('\n'), gene_name))
            continue
        if arr[gene_name_column] != '':
            output.write(line)
        else:
            arr[gene_name_column] = gene_name
            output.write('%s\n'%('\t'.join(arr)))

    output.close()

