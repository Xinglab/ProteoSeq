import re,os,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter peptide table by another peptide table.')
    parser.add_argument('-i', '--input_table', help='Input a peptide table.', type=str, required=True)
    parser.add_argument('-r', '--ref_table', help='Input a reference peptide table.', type=str)
    parser.add_argument('--retaining', help='Retaining peptide entries in the input table found in the given reference file. Set to False to discard peptides found in the reference. Default=True.', choices=['True', 'False'], default='True', type=str)
    parser.add_argument('--column_names', help='Retaining values of given columns in the reference file if the same columns are found in both inputs. Space splitted for multiple columns', nargs='+', type=str)
    parser.add_argument('-o','--output_table', help='Output peptide table', type=str)

    args = parser.parse_args()
    input_table = args.input_table
    ref_table = args.ref_table
    output_table = args.output_table
    if args.retaining == 'False':
        retaining = False
    else:
        retaining = True

    if args.column_names:
        column_names_list = args.column_names
    else:
        column_names_list = []

    peptide_seq_set = set()
    peptide_value_dict = {}
    header_dict = {}
    for line in open(ref_table, 'r'):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict) == 0:
            if line.startswith('#'):
                arr = line.rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        peptide_seq = arr[0]
        peptide_seq_set.add(peptide_seq)
        for column in column_names_list: #use values given in the reference table
            i = header_dict.get(column, -1)
            if i >= 0: 
                peptide_value_dict.setdefault(peptide_seq, {})[column] = arr[i]

    output_column_names = []
    output_line_list = []
    header_dict = {}
    for line in open(input_table, 'r'):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(1, len(arr)):
                header_dict[arr[i]] = i
                output_column_names.append(arr[i])
            for column in column_names_list:
                if column not in output_column_names:
                    output_column_names.append(column)
            continue
        peptide_seq = arr[0]
        value_list = [peptide_seq]
        for column in output_column_names:
            if column in column_names_list: #use values given in the reference table
                value = peptide_value_dict.get(peptide_seq, {}).get(column, '')
            else:
                value = arr[header_dict[column]]
            value_list.append(value)
        if peptide_seq in peptide_seq_set:
            if retaining:
                output_line_list.append('\t'.join(value_list))
        else:
            if not retaining:
                output_line_list.append('\t'.join(value_list))

    output = open(output_table, 'w')
    output.write('#%s\t%s\n'%('peptide', '\t'.join(output_column_names)))
    for output_line in output_line_list:
        output.write(output_line + '\n')
    output.close()

