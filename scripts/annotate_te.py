import os,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate types for the TE elements.')
    parser.add_argument('-i', '--input_file', help='Index table. Required', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output table. Required', type=str, required=True)
    parser.add_argument('-r', '--reference', help='Reference bed file of TE. Required', type=str, required=True)

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    reference = args.reference

    te_type_dict = {}
    for line in open(reference, 'r'):
        arr = line.rstrip('\n').split('\t')
        coord = (arr[0], arr[5], arr[1], arr[2])
        te_type_dict[coord] = (arr[7], arr[6])
        
    output = open(output_file, 'w')
    for line in open(input_file, 'r'):
        arr = line.rstrip('\n').split('\t')
        te_coord = (arr[6], arr[11], arr[7], arr[8])
        te_type = te_type_dict[te_coord]
        output.write('%s\t%s\n'%(line.rstrip('\n'), '\t'.join(te_type)))
        
    output.close()


