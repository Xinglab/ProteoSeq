import os,sys,argparse

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
    parser = argparse.ArgumentParser(description='Add new tags to all lines in a GTF file.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file. Attributes in the GTF will be kept in the output GTF. Required', type=str)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file. Required', type=str, required=True)
    parser.add_argument('--tag', help='Add tags to all lines in a format like \"attribute:value\", space splitted if multiple are given, default=None', nargs='+', type=str)
    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_gtf = args.output_gtf
    tag_list = args.tag

    add_d = {}
    for tag in tag_list:
        k,v = tag.split(':')
        add_d.setdefault(k, []).append(v)

    output = open(output_gtf, 'w')
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        d = merge_attributes(d, add_d)
        attribute = deparse_attributes(d)
        new_line = '%s\t%s\n'%('\t'.join(arr[:8]), attribute)
        output.write(new_line)

    output.close()

