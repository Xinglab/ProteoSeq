"""
Script Name: calculate_cpm.py
Description: Calculate CPM from the raw read count table.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import sys, argparse
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate CPM from the raw read count table.')
    parser.add_argument('-i', '--input_file', help='Input file', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=str, required=True)
    parser.add_argument('--skipped_columns', help='Skipped columns, comma splitted. default=1', default='1', type=str)

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    skipped_columns = [int(_)-1 for _ in args.skipped_columns.split(',')]

    total_reads_dict = {}
    header_dict = {}
    kept_columns = []
    sample_list = []
    for line in open(input_file, 'r'):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict) == 0:
            for i in range(len(arr)):
                header_dict[arr[i]] = i
                if i in skipped_columns:
                    kept_columns.append(arr[i])
                else:
                    sample_list.append(arr[i])
            continue
        for sample in sample_list:
            exp = float(arr[header_dict[sample]])
            total_reads_dict[sample] = total_reads_dict.get(sample, 0) + exp

    print('%s\t%s'%('Sample', 'Total reads'))
    for sample in sample_list:
        total_reads = total_reads_dict[sample]
        print('%s\t%d'%(sample, total_reads))

    output = open(output_file, 'w')
    header_dict = {}
    for line in open(input_file,'r'):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict) == 0:
            output.write(line)
            for i in range(len(arr)):
                header_dict[i] = arr[i]
            continue
        value_list = []
        for i in range(len(arr)):
            if i in skipped_columns:
                value_list.append(arr[i])
                continue
            sample = header_dict[i]
            total_reads = total_reads_dict[sample]
            if total_reads == 0:
                cpm = 0
            else:
                cpm = float(arr[i]) * 1000000 / total_reads
            value_list.append('%.2f'%cpm)
        n = output.write('%s\n'%('\t'.join(value_list)))

    output.close()

