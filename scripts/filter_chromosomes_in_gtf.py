import os,sys,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select transcripts encoded in certain chromosomes in a given GTF.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file', type=str, required=True)
    parser.add_argument('--output_discarded', help='Output additional GTF file of discarded lines', type=str) 
    parser.add_argument('--chromosome', help='List of chromosomes kept, comma splitted. Set off to include only primary chromosomes (1-22, X&Y)', type=str)

    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_gtf = args.output_gtf
    if args.output_discarded:
        output_discarded = args.output_discarded
    else:
        output_discarded = None

    if args.chromosome:
        chromosome_list = args.chromosome.split(',')
    else:
        chromosome_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]


    output = open(output_gtf, 'w')
    if output_discarded:
        output2 = open(output_discarded, 'w')
    else:
        output2 = None

    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[0] in chromosome_list:
            output.write(line)
        else:
            if output2:
                output2.write(line)
            else:
                print('Exclude line: %s'%(line.rstrip('\n')))

    output.close()
    if output2:
        output2.close()

