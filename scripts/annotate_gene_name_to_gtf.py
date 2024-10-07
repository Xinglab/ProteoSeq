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
            string.append('%s \"%s\"'%(k,v))
    return '; '.join(string)+';'

def merge_attributes(d1, d2):
    # keep the attribute value in d1 if the attribute exists in both dicts
    d = d1.copy()
    for k,v in d2.items():
        if not d.get(k, None):
            d[k] = v
    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate gene name by the gene id to a GTF by those in another GTF.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file to annotate. Required', type=str, required=True)
    parser.add_argument('-a', '--annotation_gtf', help='Annotation GTF with features. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF. Required', type=str, required=True)
    parser.add_argument('-v', '--verbose', help='Set to zero to silent stdout skipped lines, default=0', type=int, choices=[0, 1], default=0)

    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    annotation_gtf = args.annotation_gtf
    output_gtf = args.output_gtf
    verbose = args.verbose
    
    gene_name_dict = {}
    for line in open(annotation_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] != 'gene':
            continue
        d = parse_attributes(arr[8])
        gene_id = d.get('gene_id', [''])[0]
        if gene_id == '':
            continue
        gene_name = d.get('gene_name', [''])[0]
        if gene_name == '':
            continue
        gene_name_dict.setdefault(gene_id, {}).setdefault('gene_name', set()).add(gene_name)

    output = open(output_gtf, 'w')
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        d = parse_attributes(arr[8])
        gene_id = d.get('gene_id', [''])[0]
        if gene_id == '':
            if verbose == 1:
                print('Excluded line in annotation GTF (no gene id): %s'%(line.rstrip('\n')))
            continue
        add_d = gene_name_dict.get(gene_id, {})
        d = merge_attributes(d, add_d)
        attribute = deparse_attributes(d)
        new_line = '%s\t%s\n'%('\t'.join(arr[:8]), attribute)
        output.write(new_line)

    output.close()



