import os,argparse,re
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter significant PSMs in PIN file by output of Percolator.')
    parser.add_argument('-i', '--input_file', help='Input files with PSMs. Required', type=str, required=True)
    parser.add_argument('-p', '--input_pin', help='Input PIN files with MS2 spectra intensity. Required', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output PIN file', type=str)
    args = parser.parse_args()
    
    input_file = args.input_file
    input_pin = args.input_pin
    output_file = args.output_file

    sig_spectra_set = set()
    for line in open(input_file, 'r'):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('PSMId'):
            continue
        sig_spectra_set.add(arr[0])

    output = open(output_file, 'w')
    for line in open(input_pin, 'r'):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('DefaultDirection'):
            output.write(line)
            continue
        if line.startswith('SpecId'):
            output.write(line)
            continue
        specid = arr[0] #.split('_')[2]
        if specid in sig_spectra_set:
            output.write(line)

    output.close()


