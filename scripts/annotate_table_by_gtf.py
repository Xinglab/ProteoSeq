"""
Script Name: annotate_table_by_gtf.py
Description: Annotate index table with information from GTF.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,sys,argparse

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
    parser = argparse.ArgumentParser(description='Annotate index table with information from GTF.')
    parser.add_argument('-i', '--input_table', help='Input table protein index to original ids. Required', type=str, required=True)
    parser.add_argument('-g', '--annotation_gtf', help='Input GTF for annotation. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table with annotation retrieved from GTF. Required', type=str, required=True)
    args = parser.parse_args()

    input_table = args.input_table
    annotation_gtf = args.annotation_gtf
    output_table = args.output_table

    attribute_list = ['ORF_type', 'transcript_id', 'transcript_name', 'gene_id', 'gene_name']

    tx_attributes_dict = {}
    for line in open(annotation_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                #print('Excluded line in annotation GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            for key in attribute_list:
                value = d.get(key, [''])[0]
                tx_attributes_dict.setdefault(tx_id, {})[key] = value

    new_header_list = []
    header_dict = {}
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(len(arr)):
                new_header_list.append(arr[i])
                header_dict[arr[i]] = i
        else:
            break

    for key in attribute_list:
        if key not in new_header_list:
            new_header_list.append(key)

    output = open(output_table, 'w')
    header_rev_dict = {}
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(len(arr)):
                header_rev_dict[i] = arr[i]
            output.write('#%s\n'%('\t'.join(new_header_list)))
            continue
        arr = line.rstrip('\n').split('\t')
        if header_dict.get('transcript_id', -1) >= 0:
            tx_id = arr[header_dict['transcript_id']]
        else:
            tx_id = ''
        if tx_id == '':
            tx_id = arr[header_dict['seq_id']]
        new_value_list = []
        for i in range(len(arr)):
            if arr[i] != '':
                new_value_list.append(arr[i])
            else:
                header = header_rev_dict[i]
                if header not in attribute_list:
                    new_value_list.append(arr[i])
                else:
                    new_value = tx_attributes_dict.get(tx_id, {}).get(header, '')
                    new_value_list.append(new_value)
        for i in range(len(arr), len(new_header_list)):
            key = new_header_list[i]
            new_value = tx_attributes_dict.get(tx_id, {}).get(key, '')
            new_value_list.append(new_value)
        output.write('%s\n'%('\t'.join(new_value_list)))

    output.close()

