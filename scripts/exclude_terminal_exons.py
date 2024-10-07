import os,sys,argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Exclude terminal exons in a BED file.')
    parser.add_argument('-i', '--input_bed', help='Input BED file. Required', type=str)
    parser.add_argument('-o', '--output_bed', help='Output BED file. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_bed = args.input_bed
    output_bed = args.output_bed

    d = {}
    tx_set = set()
    tx_list = []
    for line in open(input_bed, 'r'):
        arr = line.rstrip('\n').split('\t')
        if arr[3] not in tx_set:
            tx_set.add(arr[3])
            tx_list.append(arr[3])
        d.setdefault(arr[3], []).append(line)
        
    output = open(output_bed, 'w')
    for tx_id in tx_list:
        line_list = d[tx_id]
        for line in line_list[1:-1]:
            n = output.write(line)
    output.close()


