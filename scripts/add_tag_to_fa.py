import os,sys,argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add new tags to ID lines in a FASTA file.')
    parser.add_argument('-i', '--input_fa', help='Input FASTA file. Required', type=str)
    parser.add_argument('-o', '--output_fa', help='Output FASTA file. Required', type=str, required=True)
    parser.add_argument('--tag', help='Add tags to the ID line in a FASTA file in a format like \"attribute:value\", space splitted if multiple are given, default=None', nargs='+', type=str)
    args = parser.parse_args()
    
    input_fa = args.input_fa
    output_fa = args.output_fa
    tag_list = args.tag

    add_string = ' '.join(tag_list)

    output = open(output_fa, 'w')
    for line in open(input_fa, 'r'):
        if line.startswith('>'):
            new_line = '%s %s\n'%(line.rstrip('\n'), add_string)
            output.write(new_line)
        else:
            output.write(line)

    output.close()

