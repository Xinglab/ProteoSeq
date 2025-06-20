"""
Script Name: convert_expression_to_log2.py
Description: Convert expression to log2.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,sys,argparse
import math

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert expression to log2.')
    parser.add_argument('-i', '--input_file', help='Input tabular file with expression. Required', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file. Required', type=str, required=True)
    parser.add_argument('--columns', help='Indicate columns excluded from calculation. Comma splitted. Default=1,2', type=str, default='1,2')
    parser.add_argument('--base', help='Log base. Default=2.', type=str, default='2')
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    column_excluded = [int(_)-1 for _ in args.columns.split(',')]
    if args.base == 'e':
        base = math.e
    else:
        base = int(args.base)

    output = open(output_file, 'w')
    for line in open(input_file, 'r'):
        if line.startswith('#'):
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        out_list = []
        for column in range(len(arr)):
            if column in column_excluded:
                out_list.append(arr[column])
            else:
                exp = float(arr[column])
                if exp==0:
                    out_list.append('0')
                else:
                    log = math.log(exp+1, base) #Add one
                    out_list.append('%.2f'%log)
        output.write('\t'.join(out_list)+'\n')
    output.close()

