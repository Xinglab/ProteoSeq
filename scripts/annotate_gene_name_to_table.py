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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate index table with gene name from GTF.')
    parser.add_argument('-i', '--input_table', help='Input table protein index to original ids. Required', type=str, required=True)
    parser.add_argument('-g', '--annotation_gtf', help='Input GTF for annotation. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table with annotation retrieved from GTF. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_table = args.input_table
    annotation_gtf = args.annotation_gtf
    output_table = args.output_table

    gene_name_dict = {}
    for line in open(annotation_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                #print('Excluded line in annotation GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            gene_id_list = d.get('gene_id', [''])
            gene_name_list = d.get('gene_name', [''])
            for gene_id in gene_id_list:
                if gene_id == '':
                    continue
                gene_name_dict[gene_id] = gene_name_dict.setdefault(gene_id, []) + gene_name_list
    
    gene_id_idx = -1
    gene_name_idx = -1
    header_dict = {}
    output = open(output_table, 'w')
    for line in open(input_table, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n').split('\t')
            for i in range(len(arr)):
                if arr[i] == 'gene_id':
                    gene_id_idx = i
                elif arr[i] == 'gene_name':
                    gene_name_idx = i
            if gene_id_idx < 0:
                exit('No gene id is found in the input table.')
            if gene_name_idx < 0:
                output.write('%s\t%s\n'%(line.rstrip('\n'), 'gene_name'))
            else:
                output.write(line)
            continue
        arr = line.rstrip('\n').split('\t')
        if gene_name_idx >= 0 and arr[gene_name_idx] != '':
            output.write(line)
            continue
        gene_id = arr[gene_id_idx]
        gene_name_list = gene_name_dict.get(gene_id, [''])
        gene_name = ','.join(set(gene_name_list))
        if gene_name_idx >= 0:
            arr[gene_name_idx] = gene_name
            output.write('%s\n'%('\t'.join(arr)))
        else:
            output.write('%s\t%s\n'%(line.rstrip('\n'), gene_name))

    output.close()

