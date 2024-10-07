import os,sys,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse ID table')
    parser.add_argument('-i', '--input_table', help='Input table. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table. Required', type=str, required=True)
    args = parser.parse_args()

    input_table = args.input_table
    output_table = args.output_table

    output = open(output_table, 'w')
    n = 0
    for line in open(input_table, 'r'):
        if n == 0:
            n += 1
            _ = output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        new_arr = []
        for i in range(1, len(arr)):
            arr_list = arr[i].split(';')
            if len(set(arr_list)) == 1:
                new_arr.append(arr_list[0])
            else:
                new_arr.append(arr[i])
        _ = output.write('\t'.join(new_arr) + '\n')

    output.close()
