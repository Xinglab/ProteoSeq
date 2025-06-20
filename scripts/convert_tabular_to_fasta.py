"""
Script Name: convert_tabular_to_fasta.py
Description: Convert tabular file to fasta file of sequences. Sequences will be reindexed from 1-N.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,sys,argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert tabular file to fasta file of sequences. Sequences will be reindexed from 1-N.')
    parser.add_argument('-i', '--input_file', help='Input tabular file with sequences. Required', type=str)
    parser.add_argument('-o', '--output_file', help='Output fasta file. Required', type=str, required=True)
    parser.add_argument('--column', help='Indicate which column includes sequences in the input. default=1', type=int, default=1)
    parser.add_argument('--seq_prefix', help='Add prefix to the index of each sequence, the name of the sequence will be >Seq_prefix_Index. default=None', type=str)
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    column = int(args.column)

    if args.seq_prefix:
        seq_prefix = args.seq_prefix
    else:
        seq_prefix = None

    seq_idx_dict = {}
    idx = 1
    for line in open(input_file, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        seq = arr[column-1]
        if seq_idx_dict.get(seq, 0) == 0:
            seq_idx_dict[seq] = idx
            idx += 1

        
    output = open(output_file, 'w')
    for seq, idx in seq_idx_dict.items():
        if seq_prefix:
            seq_id = '%s_%d'%(seq_prefix, idx)
        else:
            seq_id = '%d'%(idx)
        output.write('>%s\n%s\n'%(seq_id, seq))

    output.close()

