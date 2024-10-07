import os, argparse, re

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
    parser = argparse.ArgumentParser(description='Select lines with required tags in a GTF file.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file. Attributes in the GTF will be kept in the output GTF. Required', type=str)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file. Required', type=str, required=True)
    parser.add_argument('--output_discarded', help='Output file of discarded lines.', type=str)
    parser.add_argument('--tag', help='Tags to be selected in a format like \"attribute:value\", space splitted if multiple are given, default=None', nargs='+', type=str)
    parser.add_argument('--keep', help='Keep a line if it contains ANY or ALL required tags', choices=['All', 'Any'], default='All', type=str)
    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_gtf = args.output_gtf
    tag_list = args.tag
    keep = args.keep
    if args.output_discarded:
        output_file_discarded = args.output_discarded
    else:
        output_file_discarded = None

    select_d = {}
    for tag in tag_list:
        k,v = tag.split(':')
        select_d.setdefault(k, []).append(v)
    #print(select_d)

    if output_file_discarded:
        output_discarded = open(output_file_discarded, 'w')
    else:
        output_discarded = None

    if not os.path.exists(input_gtf):
        exit('Input GTF not exists!')

    output = open(output_gtf, 'w')
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        if keep == 'Any':
            for k, v_required_list in select_d.items():
                for v_required in v_required_list:
                    if v_required in d.get(k, []):
                        output.write(line) #Output the line if it has ANY of required tags
                        break #Break to next line
                else:
                    continue #Continue to next attribute
                break #Having ANY of required tags, break to next line
            else: #Finished looping all required tags w/o breaking
                if output_discarded:
                    output_discarded.write(line)
        elif keep == 'All':
            for k, v_required_list in select_d.items():
                for v_required in v_required_list:
                    if v_required not in d.get(k, []):
                        if output_discarded:
                            output_discarded.write(line)
                        break #Not having all required tags
                else:
                    continue
                break #Not having all required tags
            else: #Finished looping all required tags w/o breaking
                output.write(line)
        
    output.close()

