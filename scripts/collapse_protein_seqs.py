import os,sys,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collapse identical sequence and reindex unique sequences.')
    parser.add_argument('-i', '--input_file_list', help='Input FASTA files, space splitted. Required', nargs='+', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output FASTA. Required', type=str, required=True)
    parser.add_argument('-ot', '--output_table', help='Output table of protein index to original ids. Required', type=str, required=True)
    parser.add_argument('--add_tags', help='Add tag in a format like \"attribute:tag1,tag2,tag3\" as a column to each of the input file. If not empty, number tags (comma splitted) must be same as the number of input_file_list, default=None', type=str)
    args = parser.parse_args()
    
    input_file_list = args.input_file_list
    output_file = args.output_file
    output_table_file = args.output_table
    if args.add_tags:
        arr = args.add_tags.split(':')
        if len(arr) != 2:
            exit()
        attribute = arr[0]
        add_tag_list = arr[1].split(',')
        if len(add_tag_list) != len(input_file_list):
            exit('Number of add_tags is different from the number of input files!')
    else:
        attribute = None
        add_tag_list = None
    
    protein_seq_dict_all_files = {}
    file_tag_dict = {}
    idx = 0
    for idx in range(len(input_file_list)):
        input_file = input_file_list[idx]
        if add_tag_list:
            file_tag_dict[input_file] = add_tag_list[idx]
        else:
            file_tag_dict[input_file] = None
        protein_seq_dict = {}
        for line in open(input_file, 'r'):
            if line[0] == '>':
                seq_id = line.rstrip('\n')[1:]
                if seq_id.startswith('ENST'): #Ensembl format transcript ID
                    seq_id = seq_id.split()[0]
                elif seq_id.startswith('ESPRESSO'): #ESPRESSO format transcript ID
                    seq_id = seq_id.split()[0]
                elif (seq_id.startswith('ENSP')) and (seq_id.find('|') >= 0): #GENCODE pc_translation format transcript ID
                    seq_id = seq_id.split('|')[1]
                elif seq_id.find('|') >= 0: #UniProt format protein ID
                    seq_id = seq_id.split('|')[1]
                else: #Other customized inputs, use full header
                    #seq_id = seq_id.split()[0]
                    #print(seq_id)
                    pass
                #if protein_seq_dict.get(seq_id, None):
                #    print('ID %s duplicates in input %s'%(seq_id, input_file))
            else:
                protein_seq_dict[seq_id] = protein_seq_dict.get(seq_id, '') + line.rstrip('\n').upper()
        protein_seq_dict_all_files[input_file] = protein_seq_dict.copy()

    protein_seq_dict = {}
    seq_tag_dict = {}
    for input_file, d in protein_seq_dict_all_files.items():
        add_tag = file_tag_dict[input_file]
        for seq_id, seq in d.items():
            seq_tag_dict.setdefault(seq_id, set()).add(add_tag)
            seq2 = protein_seq_dict.get(seq_id, None)
            if seq2 != None and seq != seq2:
                print('ID %s has different sequences in references:\nseq1: %s\nseq2: %s'%(seq_id, seq2, seq1))
            protein_seq_dict[seq_id] = seq

    protein_seq_all_tx_dict = {}
    for seq_id, protein_seq in protein_seq_dict.items():
        protein_seq_all_tx_dict.setdefault(protein_seq, []).append(seq_id)
    print('Total protein count: %d'%len(protein_seq_all_tx_dict))

    output = open(output_file, 'w')
    output_table = open(output_table_file, 'w')
    if attribute != None:
        output_table.write('#protein_idx\tseq_id\t%s\n'%(attribute))
    else:
        output_table.write('#protein_idx\tseq_id\n')
    i = 0
    for protein_seq in sorted(protein_seq_all_tx_dict, key=lambda k:len(protein_seq_all_tx_dict[k]), reverse=True):
        i += 1
        seq_id_list = protein_seq_all_tx_dict[protein_seq]
        output.write('>Protein_%d\n%s\n'%(i, protein_seq))
        for seq_id in seq_id_list:
            if attribute != None:
                add_tag = ','.join(seq_tag_dict.get(seq_id, set('')))
                output_table.write('Protein_%d\t%s\t%s\n'%(i, seq_id, add_tag))
            else:
                output_table.write('Protein_%d\t%s\n'%(i, seq_id))

    output.close()
    output_table.close()
