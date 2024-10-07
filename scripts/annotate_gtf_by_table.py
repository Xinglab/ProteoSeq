import os,sys,argparse

def parse_attributes(string):
    d = {}
    for i in string.split(';'):
        if i == '':
            continue
        i_list = i.strip().split(' ')
        if len(i_list) < 2:
            print('Wrong attributes: %s'%i)
            continue
        d.setdefault(i_list[0], []).append(i_list[1].strip('"'))
    return d


def deparse_attributes(d):
    string = []
    for k,v_list in d.items():
        for v in v_list:
            string.append('%s \"%s\"'%(k, v))
    return '; '.join(string)+';'


def merge_attributes(d1, d2):
    # keep the attribute value in d1 if the attribute exists in both dicts
    d = d1.copy()
    for k,v in d2.items():
        if not d.get(k, None):
            d[k] = v
    return d


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate additional information from text table to GTF.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF. Required', type=str, required=True)
    parser.add_argument('-t', '--input_table', help='Input table with additional information. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF. Required', type=str, required=True)
    parser.add_argument('--columns', help='Columns in input_table containing information to be added to GTF, comma splitted. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    input_table = args.input_table
    output_gtf = args.output_gtf
    columns = [int(_) for _ in args.columns.split(',')]

    new_attribute_dict = {}
    header_dict = {}
    for line in open(input_table, 'r'):
        arr = line.rstrip('\n').split('\t')
        if len(header_dict) == 0:
            for i in range(len(arr)):
                header_dict[i] = arr[i].lower()
            continue
        tx_id = arr[0]
        for c in columns:
            key = header_dict[c-1]
            values = arr[c-1].split(',')
            for value in values:
                if (value != 'NA') and (value != ''):
                    new_attribute_dict.setdefault(tx_id, {}).setdefault(key, set()).add(value)

    output = open(output_gtf, 'w')
    for line in open(input_gtf,'r'):
        if line[0] == '#':
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        tx_id = d.get('transcript_id', [''])[0]
        if tx_id == '':
            output.write(line)
            continue
        d_bed = new_attribute_dict.get(tx_id, {})
        d_merged = merge_attributes(d, d_bed)
        attribute = deparse_attributes(d_merged)
        output.write('%s\t%s\n'%('\t'.join(arr[:8]), attribute))

    output.close()


