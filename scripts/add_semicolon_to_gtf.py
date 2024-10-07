import os,sys,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add missing semicolon to the end of lines in a GTF file.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file. Required', type=str)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file. Required', type=str, required=True)

    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_gtf = args.output_gtf

    output = open(output_gtf, 'w')
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip().split('\t')
        if not arr[8].endswith(';'):
            output.write('%s\t%s;\n'%('\t'.join(arr[:8]), arr[8]))
        else:
            output.write(line)

    output.close()

