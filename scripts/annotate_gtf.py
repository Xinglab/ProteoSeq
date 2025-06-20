"""
Script Name: annotate_gtf.py
Description: Annotate features and attributes to a GTF by those in another GTF.
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
            continue
        d.setdefault(i_list[0], []).append(i_list[1].strip('"'))
    return d


def deparse_attributes(d):
    string = []
    for k,v_list in d.items():
        for v in v_list:
            string.append('%s \"%s\"'%(k,v))
    return '; '.join(string)


def merge_attributes(d1, d2):
    # keep the attribute value in d1 if the attribute exists in both dicts
    d = d1.copy()
    for k,v in d2.items():
        if not d.get(k, None):
            d[k] = v
    return d


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate features and attributes to a GTF by those in another GTF.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file to annotate. Required', type=str, required=True)
    parser.add_argument('-a', '--annotation_gtf', help='Annotation GTF with features. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF. Required', type=str, required=True)
    parser.add_argument('--selected_features', help='Features in the annotation GTF to be added to the output, comma splitted, default=CDS', type=str, default='CDS')
    parser.add_argument('-v', '--verbose', help='Set to zero to silent stdout skipped lines, default=0', type=int, choices=[0, 1], default=0)

    args = parser.parse_args()
    input_gtf = args.input_gtf
    annotation_gtf = args.annotation_gtf
    output_gtf = args.output_gtf
    selected_features = args.selected_features.split(',')
    verbose = args.verbose

    tx_attributes_dict = {}
    tx_features_dict = {}
    for line in open(annotation_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] not in ['transcript'] + selected_features:
            continue
        d = parse_attributes(arr[8])
        tx_id = d.get('transcript_id', [''])[0]
        if tx_id == '':
            if verbose == 1:
                print('Excluded line in annotation GTF (no transcript id): %s'%(line.rstrip('\n')))
            continue
        if arr[2] == 'transcript':
            tx_attributes_dict[tx_id] = d
        elif arr[2] in selected_features:
            tx_features_dict.setdefault(tx_id, []).append(line)

    tx_annotated_lines = {}
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] in selected_features:
            continue #Old feature lines in the input file will be discarded
        d1 = parse_attributes(arr[8])
        tx_id = d1.get('transcript_id', [''])[0]
        if tx_id == '':
            if verbose == 1:
                print('Excluded line in input GTF (no transcript id): %s'%(line.rstrip('\n')))
            continue
        d2 = tx_attributes_dict.get(tx_id, None)
        if not d2:
            tx_annotated_lines.setdefault(tx_id, []).append(line)
            continue
        d = merge_attributes(d1, d2)
        attribute = deparse_attributes(d)
        new_line = '%s\t%s\n'%('\t'.join(arr[:8]), attribute)
        tx_annotated_lines.setdefault(tx_id, []).append(new_line)


    output = open(output_gtf, 'w')
    for tx_id, new_line_list in tx_annotated_lines.items():
        for new_line in new_line_list:
            output.write(new_line)
        feature_list = tx_features_dict.get(tx_id, [])
        if len(feature_list) > 0:
            for feature in feature_list:
                output.write(feature)

    output.close()


