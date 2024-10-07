import os,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select significant PSMs in Percolator output.')
    parser.add_argument('-i', '--input_file', help='Input file from Percolator. Required', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=str)
    parser.add_argument('-f', '--fdr', help='Fdr threshold, default=0.05', type=float, default='0.05')
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    fdr_thred = float(args.fdr)

    output = open(output_file, 'w')
    header_dict = {}
    for line in open(input_file, 'r'):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict)==0:
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            output.write(line)
            continue
        qvalue = float(arr[header_dict['q-value']])
        if qvalue >= fdr_thred:
            continue
        output.write(line)

    output.close()



