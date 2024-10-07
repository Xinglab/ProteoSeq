import re,os,argparse
import numpy as np
    
def get_protein_id_to_seq_dict(input_file):
    protein_id_to_seq_dict = {}
    for line in open(input_file, 'r'):
        if line[0]=='>':
            protein_id = line.rstrip('\n')[1:]
        else:
            protein_id_to_seq_dict[protein_id] = protein_id_to_seq_dict.get(protein_id, '') + line.rstrip('\n').upper()
    return protein_id_to_seq_dict


def get_peptide_mapping_dicts(input_file, exclude_list=None):
    peptide_count_dict = {}
    sample_list = []
    header_dict = {}
    abun_idx = 1
    for line in open(input_file, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
                if arr[i].startswith('total_'):
                    if abun_idx > 1:
                        continue
                    abun_idx = i
            sample_list = arr[abun_idx:]
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[header_dict['peptide']]
        peptide_count_dict[peptide_seq] = arr[abun_idx:]
    return peptide_count_dict,sample_list


def outf_table(peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_ends, output_file, append=False):
    if not append:
        outf = open(output_file, 'w')
        out_line = '#%s\t%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id', '\t'.join(sample_list))
        outf.write(out_line)
    else:
        outf = open(output_file, 'a')
    n = 0
    for peptide_seq in sorted(peptide_all_mapping_dict, key=lambda k:peptide_count_dict[k], reverse=True):
        tx_list = peptide_all_mapping_dict[peptide_seq]
        count_list = peptide_count_dict[peptide_seq]
        start_list = peptide_all_mapping_dict_starts[peptide_seq]
        end_list = peptide_all_mapping_dict_ends[peptide_seq]
        out_line = '%s\t%s\t%s\t%s\t%s\n'%(peptide_seq, ','.join(start_list), ','.join(end_list), ','.join(tx_list), '\t'.join(count_list))
        outf.write(out_line)
    outf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter peptides in exsiting FASTA.')
    parser.add_argument('-pep', '--peptide_file', help='Input peptide table. Required.', required=True)
    parser.add_argument('-ref', '--ref_file', help='Input FASTA file of reference protein sequences. Required', required=True)
    parser.add_argument('-o2', '--output_unmapped', help='Output peptides not in reference FASTA', type=str)
    parser.add_argument('-o1', '--output_mapped', help='Output peptides in reference FASTA, mappings are added as additional columns.', type=str)
    args = parser.parse_args()
    
    ref_file = args.ref_file
    peptide_file = args.peptide_file
    if args.output_unmapped:
        output_file_unmapped = args.output_unmapped
    else:
        output_file_unmapped = None
        
    if args.output_mapped:
        output_file_mapped = args.output_mapped
    else:
        output_file_mapped = None

    protein_id_to_seq_dict = get_protein_id_to_seq_dict(ref_file)
    protein_seq_all_id_dict = {}
    for seq_id, protein_seq in protein_id_to_seq_dict.items():
        protein_seq_all_id_dict.setdefault(protein_seq, []).append(seq_id)

    #print('Entries in reference: %d'%(len(protein_id_to_seq_dict)))
    #print('Unique proteins in reference: %d'%(len(protein_seq_all_id_dict)))

    peptide_count_dict,sample_list = get_peptide_mapping_dicts(peptide_file) #, exclude_list=set(novel_peptide_count_dict.keys()))

    peptide_all_mapping_dict = {}
    peptide_all_mapping_dict_starts = {}
    peptide_all_mapping_dict_ends = {}
    for peptide_seq in peptide_count_dict.keys():
        all_tx_list = []
        all_tx_start_list = []
        all_tx_end_list = []
        for protein_seq, all_seq_list in protein_seq_all_id_dict.items():
            for match in re.finditer(peptide_seq, protein_seq):
                start = match.start()+1 #Python mode
                end = match.end()
                for seq_id in all_seq_list: 
                    all_tx_list.append(seq_id)
                    all_tx_start_list.append(str(start))
                    all_tx_end_list.append(str(end))
                    #break
                #break
        if len(all_tx_list) > 0:
            peptide_all_mapping_dict[peptide_seq] = all_tx_list
            peptide_all_mapping_dict_starts[peptide_seq] = all_tx_start_list
            peptide_all_mapping_dict_ends[peptide_seq] = all_tx_end_list

    unmapped_peptide_set = set(peptide_count_dict) - set(peptide_all_mapping_dict)

    print('%s\t%s\t%s'%('#sample', 'mapped_peptides', 'unmapped_peptides'))
    print('%s\t%d\t%d'%(peptide_file, len(peptide_all_mapping_dict), len(unmapped_peptide_set)))

    if output_file_mapped:
        outf_table(peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_ends, output_file_mapped)

    if output_file_unmapped:
        outf = open(output_file_unmapped, 'w')
        for line in open(peptide_file, 'r'):
            if line.startswith('#'):
                outf.write(line)
                continue
            arr = line.rstrip('\n').split('\t')
            if arr[0] in unmapped_peptide_set:
                outf.write(line)
        outf.close()

