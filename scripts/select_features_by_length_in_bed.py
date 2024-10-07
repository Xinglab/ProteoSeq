import os,sys,argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select features with certain length in a BED file.')
    parser.add_argument('-i', '--input_bed', help='Input BED file. Required', type=str, required=True)
    parser.add_argument('-o', '--output_bed', help='Output BED file. Required', type=str, required=True)
    parser.add_argument('-min', '--min_length', help='Minimum length (bp) of the feature. Set to zero to skip this filter. Default=0.', type=int, default=0)
    parser.add_argument('-max', '--max_length', help='Maximum length (bp) of the feature. Set to zero to skip this filter. Default=0.', type=str, default=0)
    parser.add_argument('-s', '--start_column', help='Column with start coordinate in the BED for filtering. Set to zero to skip this filter. Default=2', type=int, default=0)
    parser.add_argument('-e', '--end_column', help='Column with end coordinate in the BED for filtering. Set to zero to skip this filter. Default=3', type=int, default=0)
    parser.add_argument('-l', '--overlap_column', help='Column with overlapped length in the BED for filtering. Set to zero to skip this filter. Default=0', type=int, default=0)
    parser.add_argument('-v', '--verbose', help='Set to zero to silent output on screen. Default=0', type=int, choices=[0,1], default=0)
    args = parser.parse_args()
    
    input_bed = args.input_bed
    output_bed = args.output_bed
    min_length = int(args.min_length)
    max_length = int(args.max_length)
    start_column = int(args.start_column) - 1
    end_column = int(args.end_column) - 1
    overlap_column = int(args.overlap_column) - 1
    verbose = int(args.verbose)

    output = open(output_bed, 'w')
    for line in open(input_bed, 'r'):
        arr = line.rstrip('\n').split('\t')
        if start_column > -1:
            start = int(arr[start_column]) + 1 #BED format
        else:
            start = 0
        if end_column > -1:
            end = int(arr[end_column])
        else:
            end = 0
        if (start > 0) and (end > 0):
            length = end - start + 1
            if (min_length > 0) and (length < min_length):
                if verbose == 1:
                    print('Excluded:%s'%line)
                continue
            if (max_length > 0) and (length > max_length):
                if verbose == 1:
                    print('Excluded:%s'%line)
                continue
        if overlap_column > -1:
            overlap = int(arr[overlap_column])
        else:
            overlap = 0
        if overlap > 0:
            if (min_length > 0) and (overlap < min_length):
                if verbose == 1:
                    print('Excluded:%s'%line)
                continue
            if (max_length > 0) and (overlap > max_length):
                if verbose == 1:
                    print('Excluded:%s'%line)
                continue
        output.write(line)

    output.close()

