import re,os,argparse

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
    parser = argparse.ArgumentParser(description='Filter peptide table by another peptide table.')
    parser.add_argument('-i', '--input_table', help='Input a peptide table.', type=str, required=True)
    parser.add_argument('-g', '--ref_gtf', help='Input a reference GTF with peptides (peptide_seq must be in the tags).', type=str, required=True)
    parser.add_argument('--retaining', help='Retaining peptide entries in the input table found in the given reference file. Set to False to discard peptides found in the reference. Default=True.', choices=['True', 'False'], default='True', type=str)
    parser.add_argument('--column_names', help='Retaining values of given columns in the reference file if the same columns are found in both inputs. Space splitted for multiple columns', nargs='+', type=str)
    parser.add_argument('-o','--output_table', help='Output peptide table', type=str)

    args = parser.parse_args()
    input_table = args.input_table
    ref_gtf = args.ref_gtf
    output_table = args.output_table
    if args.retaining == 'False':
        retaining = False
    else:
        retaining = True

    peptide_seq_set = set()
    for line in open(ref_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        for peptide_seq in d.get('peptide_seq', []):
            peptide_seq_set.add(peptide_seq)

    if not os.path.exists(input_table):
        exit('Input table not exists!')
        
    output = open(output_table, 'w')
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[0]
        if peptide_seq in peptide_seq_set:
            if retaining:
                output.write(line)
        else:
            if not retaining:
                output.write(line)

    output.close()

