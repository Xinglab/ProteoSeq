"""
Script Name: filter_gtf_by_cpm.py
Description: Filter transcripts with their expression in given samples.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,argparse
import numpy as np

def parse_attributes(string):
    d = {}
    for i in string.split(';'):
        if i == '':
            continue
        i_list = i.strip().split(' ')
        if len(i_list) < 2:
            print(string)
            continue
        d.setdefault(i_list[0], []).append(i_list[1].strip('"'))
    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter transcripts with their expression in given samples.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file. Required', type=str, required=True)
    parser.add_argument('-e', '--exp_file', help='Tabular table of transcripts expression. Columns: samples. Rows: transcripts. Required', type=str, required=True)
    parser.add_argument('--sample', help='Add to indicate in which sample the expression should be used for filtering. Use comma to split samples if multiple samples are input. Default is all samples', type=str)
    parser.add_argument('--skipped_columns', help='Skipped columns. Comma splitted.', type=str)
    parser.add_argument('--cutoff', help='Expression cutoff, default=0', default=0, type=float)
    parser.add_argument('--value', help='When multiple samples are given to --sample, set to choose how expression level in various samples should be calculated for filtering. default=max', type=str, default='max', choices=['max', 'min', 'avg', 'sum', 'median'])
    parser.add_argument('--include_equal', help='Add this flag to include equal to in addition to greater than, default=True', type=str, default='True', choices=['True', 'False'])
    args = parser.parse_args()

    input_gtf = args.input_gtf
    output_gtf = args.output_gtf
    exp_file = args.exp_file

    cutoff = float(args.cutoff)
    if not args.sample:
        used_sample_list = []
    else:
        used_sample_list = args.sample.split(',')

    value = args.value

    if not args.skipped_columns:
        skipped_columns = []
    else:
        skipped_columns = [int(_)-1 for _ in args.skipped_columns.split(',')]

    if args.include_equal=='True':
        include_equal = True
    else:
        include_equal = False

    header_dict = {}
    tx_exp_dict = {}
    for line in open(exp_file, 'r'):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict) == 0:
            for i in range(1, len(arr)):
                if i not in skipped_columns:
                    header_dict[arr[i]] = i
            continue
        tx_id = arr[0]
        for sample, i in header_dict.items():
            exp = float(arr[i])
            if (len(used_sample_list) > 0) and (sample not in used_sample_list):
                continue
            tx_exp_dict[tx_id][sample] = tx_exp_dict.setdefault(tx_id, {}).get(sample, 0) + exp

    tx_filtered_set = set()
    for tx_id, exp_d in tx_exp_dict.items():
        exp_list = list(exp_d.values())
        if len(exp_list) == 0:
            exp = 0
        elif value == 'sum':
            exp = sum(exp_list)
        elif value == 'avg':
            exp = sum(exp_list) / len(exp_list)
        elif value == 'min':
            exp = min(exp_list)
        elif value == 'max':
            exp = max(exp_list)
        elif value == 'median':
            exp = np.median(exp_list)
        else:
            exp = sum(exp_list)
        if include_equal:
            if exp >= cutoff:
                tx_filtered_set.add(tx_id)
        else:
            if exp > cutoff:
                tx_filtered_set.add(tx_id)

    print('Input GTF: %d transcripts\nOutput GTF: %d transcripts'%(len(tx_exp_dict), len(tx_filtered_set)))
    output = open(output_gtf, 'w')
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        tx_id = d.get('transcript_id', [''])[0]
        if tx_id == '':
            continue
        if tx_id in tx_filtered_set:
            output.write(line)

    output.close()

