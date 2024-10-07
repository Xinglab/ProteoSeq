import re,os,argparse
import numpy as np
    
def get_protein_idx_to_seq_dict(input_file):
    protein_idx_to_seq_dict = {}
    for line in open(input_file, 'r'):
        if line[0]=='>':
            protein_idx = int(line.rstrip('\n')[1:].replace('Protein_', ''))
        else:
            protein_idx_to_seq_dict[protein_idx] = protein_idx_to_seq_dict.get(protein_idx, '') + line.rstrip('\n').upper()
    return protein_idx_to_seq_dict


def get_protein_idx_to_uniprot_dict(input_file):
    protein_idx_to_all_tx_dict = {}
    header_dict = {}
    for line in open(input_file, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[header_dict['ORF_type']] == 'UniProt':
            protein_idx = int(arr[0].replace('Protein_', ''))
            tx_id = arr[1]
            protein_idx_to_all_tx_dict.setdefault(protein_idx, []).append(tx_id)
    return protein_idx_to_all_tx_dict


def get_protein_idx_to_all_tx_dict(input_file):
    header_dict = {}
    for line in open(input_file, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            break
    c = header_dict.get('ORF_type', None)
    protein_idx_to_all_tx_dict = {}
    tx_orf_type_dict = {}
    for line in open(input_file, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = int(arr[0].replace('Protein_', ''))
        tx_id = arr[1].split(' ')[0]
        protein_idx_to_all_tx_dict.setdefault(protein_idx, []).append(tx_id)
        if c != None:
            tx_orf_type_dict[tx_id] = arr[c]
    if len(tx_orf_type_dict) == 0:
        return protein_idx_to_all_tx_dict, {}
    protein_orf_type_dict = {}
    for protein_idx, tx_list in protein_idx_to_all_tx_dict.items():
        for tx_id in tx_list:
            if tx_orf_type_dict[tx_id] == 'UniProt':
                protein_orf_type_dict[protein_idx] = 'UniProt'
                break
        else:
            for tx_id in tx_list:
                if tx_orf_type_dict[tx_id] == 'GENCODE': # or tx_id.startswith('ENST')
                    protein_orf_type_dict[protein_idx] = 'GENCODE'
                    break
            else:
                for tx_id in tx_list:
                    if tx_orf_type_dict[tx_id] == 'Canonical':
                        protein_orf_type_dict[protein_idx] = 'Canonical'
                        break
                else:
                    protein_orf_type_dict[protein_idx] = 'Novel'
    return protein_idx_to_all_tx_dict, protein_orf_type_dict


def get_separate_protein_seq_dicts(protein_orf_type_dict):
    protein_seq_dict_uniprot = {}
    protein_seq_dict_gencode = {}
    protein_seq_dict_canonical = {}
    protein_seq_dict_novel = {}
    for protein_idx, orf_type in protein_orf_type_dict.items():
        protein_seq = protein_idx_to_seq_dict[protein_idx]
        if orf_type == 'UniProt':
            protein_seq_dict_uniprot[protein_idx] = protein_seq
        elif orf_type == 'GENCODE':
            protein_seq_dict_gencode[protein_idx] = protein_seq
        elif orf_type == 'Canonical':
            protein_seq_dict_canonical[protein_idx] = protein_seq
        elif orf_type == 'Novel':
            protein_seq_dict_novel[protein_idx] = protein_seq
        else:
            print(protein_idx)
    return protein_seq_dict_uniprot, protein_seq_dict_gencode, protein_seq_dict_canonical, protein_seq_dict_novel


def first_uniprot_filter(peptide_to_protein_idx_dict, protein_orf_type_dict, protein_seq_dict_uniprot, protein_idx_to_uniprot_dict, report_all_hits=False):
    n_uniprot1 = 0
    selected_peptide_dict = {}
    peptide_all_mapping_dict = {}
    peptide_all_mapping_dict_starts = {}
    peptide_all_mapping_dict_dict_ends = {}
    for peptide_seq, protein_idx_list in peptide_to_protein_idx_dict.items():
        for protein_idx in protein_idx_list:
            orf_type = protein_orf_type_dict[protein_idx]
            if orf_type == 'UniProt':
                if not report_all_hits:
                    protein_seq = protein_seq_dict_uniprot[protein_idx]
                    start = protein_seq.find(peptide_seq) + 1
                    end = start + len(peptide_seq) - 1
                    peptide_all_mapping_dict[peptide_seq] = [str(protein_idx)]
                    peptide_all_mapping_dict_starts[peptide_seq] = [str(start)]
                    peptide_all_mapping_dict_dict_ends[peptide_seq] = [str(end)]
                else:
                    all_tx_list = []
                    all_tx_start_list = []
                    all_tx_end_list = []
                    for protein_idx, protein_seq in protein_seq_dict_uniprot.items():
                        start = protein_seq.find(peptide_seq) + 1
                        if start >= 1:
                            end = start + len(peptide_seq) - 1
                            tx_list2 = protein_idx_to_uniprot_dict.get(protein_idx, [str(protein_idx)])
                            all_tx_list += tx_list2
                            all_tx_start_list += [str(start)]*len(tx_list2)
                            all_tx_end_list += [str(end)]*len(tx_list2)
                    peptide_all_mapping_dict[peptide_seq] = all_tx_list
                    peptide_all_mapping_dict_starts[peptide_seq] = all_tx_start_list
                    peptide_all_mapping_dict_dict_ends[peptide_seq] = all_tx_end_list
                n_uniprot1 += 1
                break
        else:
            selected_peptide_dict[peptide_seq] = protein_idx_list
    #print('UniProt: %d'%(n_uniprot1))
    #print('Others: %d'%(len(selected_peptide_dict)))
    if not report_all_hits:
        return n_uniprot1, selected_peptide_dict
    return n_uniprot1, selected_peptide_dict, peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_dict_ends


def second_uniprot_filter(peptide_dict, protein_seq_dict_uniprot, protein_idx_to_uniprot_dict, report_all_hits=False):
    selected_peptide_dict = {}
    peptide_all_mapping_dict = {}
    peptide_all_mapping_dict_starts = {}
    peptide_all_mapping_dict_dict_ends = {}
    for peptide_seq, protein_idx_list in peptide_dict.items():
        all_tx_list = []
        all_tx_start_list = []
        all_tx_end_list = []
        for protein_idx, protein_seq in protein_seq_dict_uniprot.items():
            start = protein_seq.find(peptide_seq) + 1
            if start >= 1:
                end = start + len(peptide_seq) - 1
                tx_list2 = protein_idx_to_uniprot_dict.get(protein_idx, [str(protein_idx)])
                all_tx_list += tx_list2
                all_tx_start_list += [str(start)]*len(tx_list2)
                all_tx_end_list += [str(end)]*len(tx_list2)
                if not report_all_hits:
                    break
        if len(all_tx_list) == 0:
            selected_peptide_dict[peptide_seq] = protein_idx_list
        else:
            peptide_all_mapping_dict[peptide_seq] = all_tx_list
            peptide_all_mapping_dict_starts[peptide_seq] = all_tx_start_list
            peptide_all_mapping_dict_dict_ends[peptide_seq] = all_tx_end_list
    n_uniprot2 = len(peptide_all_mapping_dict)
    #print('UniProt: %d'%(n_uniprot2))
    #print('Others: %d'%(len(selected_peptide_dict)))
    if not report_all_hits:
        return n_uniprot2, selected_peptide_dict
    return n_uniprot2, selected_peptide_dict, peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_dict_ends


def third_gencode_filter(peptide_dict, protein_seq_dict_gencode, report_all_hits=False):
    selected_peptide_dict = {}
    peptide_all_mapping_dict = {}
    peptide_all_mapping_dict_starts = {}
    peptide_all_mapping_dict_dict_ends = {}
    for peptide_seq, protein_idx_list in peptide_dict.items():
        all_tx_list = []
        all_tx_start_list = []
        all_tx_end_list = []
        for protein_idx, protein_seq in protein_seq_dict_gencode.items():
            start = protein_seq.find(peptide_seq) + 1
            if protein_seq.find(peptide_seq) >= 0:
                end = start + len(peptide_seq) - 1
                tx_list2 = protein_idx_to_all_tx_dict[protein_idx]
                all_tx_list += tx_list2
                all_tx_start_list += [str(start)]*len(tx_list2)
                all_tx_end_list += [str(end)]*len(tx_list2)
                if not report_all_hits:
                    break
        if len(all_tx_list) == 0:
            selected_peptide_dict[peptide_seq] = protein_idx_list
        else:
            peptide_all_mapping_dict[peptide_seq] = all_tx_list
            peptide_all_mapping_dict_starts[peptide_seq] = all_tx_start_list
            peptide_all_mapping_dict_dict_ends[peptide_seq] = all_tx_end_list
    n_gencode = len(peptide_all_mapping_dict)
    #print('GENCODE: %d'%(n_gencode))
    #print('Others: %d'%(len(selected_peptide_dict)))
    if not report_all_hits:
        return n_gencode, selected_peptide_dict
    return n_gencode, selected_peptide_dict, peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_dict_ends


def fourth_canonical_filter(peptide_dict, protein_seq_dict_canonical, report_all_hits=False):
    selected_peptide_dict = {}
    peptide_all_mapping_dict = {}
    peptide_all_mapping_dict_starts = {}
    peptide_all_mapping_dict_dict_ends = {}
    for peptide_seq, protein_idx_list in peptide_dict.items():
        all_tx_list = []
        all_tx_start_list = []
        all_tx_end_list = []
        for protein_idx, protein_seq in protein_seq_dict_canonical.items():
            start = protein_seq.find(peptide_seq) + 1
            if protein_seq.find(peptide_seq) >= 0:
                end = start + len(peptide_seq) - 1
                tx_list2 = protein_idx_to_all_tx_dict[protein_idx]
                all_tx_list += tx_list2
                all_tx_start_list += [str(start)]*len(tx_list2)
                all_tx_end_list += [str(end)]*len(tx_list2)
                if not report_all_hits:
                    break
        if len(all_tx_list) == 0:
            selected_peptide_dict[peptide_seq] = protein_idx_list
        else:
            peptide_all_mapping_dict[peptide_seq] = all_tx_list
            peptide_all_mapping_dict_starts[peptide_seq] = all_tx_start_list
            peptide_all_mapping_dict_dict_ends[peptide_seq] = all_tx_end_list
    n_canonical = len(peptide_all_mapping_dict)
    #print('Other_canonical: %d'%(n_canonical))
    #print('Others: %d'%(len(selected_peptide_dict)))
    if not report_all_hits:
        return n_canonical, selected_peptide_dict
    return n_canonical, selected_peptide_dict, peptide_all_mapping_dict, peptide_all_mapping_dict_starts, peptide_all_mapping_dict_dict_ends


def get_all_mappings_of_novel_peptides(peptide_dict, protein_seq_dict_novel):
    selected_peptide_dict = {}
    peptide_all_mapping_dict = {}
    peptide_all_mapping_dict_starts = {}
    peptide_all_mapping_dict_dict_ends = {}
    for peptide_seq, protein_idx_list in peptide_dict.items():
        all_tx_list = []
        all_tx_start_list = []
        all_tx_end_list = []
        for protein_idx, protein_seq in protein_seq_dict_novel.items():
            start = protein_seq.find(peptide_seq) + 1
            end = start + len(peptide_seq) - 1 
            if start >= 1:
                tx_list2 = protein_idx_to_all_tx_dict[protein_idx]
                all_tx_list += tx_list2
                all_tx_start_list += [str(start)]*len(tx_list2)
                all_tx_end_list += [str(end)]*len(tx_list2)
        if len(all_tx_list) == 0:
            selected_peptide_dict[peptide_seq] = protein_idx_list
        else:
            peptide_all_mapping_dict[peptide_seq] = all_tx_list
            peptide_all_mapping_dict_starts[peptide_seq] = all_tx_start_list
            peptide_all_mapping_dict_dict_ends[peptide_seq] = all_tx_end_list
    n_novel = len(peptide_all_mapping_dict)
    #print('Novel: %d'%(n_novel))
    return n_novel,selected_peptide_dict,peptide_all_mapping_dict,peptide_all_mapping_dict_starts,peptide_all_mapping_dict_dict_ends


def get_peptide_mapping_dicts(input_file, exclude_list=None):
    peptide_to_protein_idx_dict = {}
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
        protein_idx = [int(_) for _ in arr[header_dict['protein_index']].split(';')]
        peptide_count_dict[peptide_seq] = arr[abun_idx:]
        #if len(exclude_list) > 0 and peptide_seq in exclude_list:
        #    continue
        peptide_to_protein_idx_dict[peptide_seq] = protein_idx
    return peptide_to_protein_idx_dict,peptide_count_dict,sample_list


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
    parser = argparse.ArgumentParser(description='Get novel peptides.')
    parser.add_argument('-pep', '--peptide_file', help='Input peptide table. Required.', required=True)
    parser.add_argument('-ref', '--ref_file', help='Input FASTA file of reference protein sequences. Required', required=True)
    parser.add_argument('-idx', '--index_file', help='Input file mapping indexes to reference protein sequences. Required', required=True)
    parser.add_argument('-o','--output_novel', help='Output novel peptides (not in UniProt nor gencode)', type=str)
    parser.add_argument('--output_uniprot', help='Output UniProt peptides', type=str)
    parser.add_argument('--output_gencode', help='Output GENCODE peptides (not in UniProt)', type=str)
    parser.add_argument('--output_canonical', help='Output other peptides mapped to given canonical reference', type=str)
    args = parser.parse_args()
    
    index_file = args.index_file
    ref_file = args.ref_file
    peptide_file = args.peptide_file
    if args.output_novel:
        output_file_novel = args.output_novel
    else:
        output_file_novel = None

    if args.output_uniprot:
        output_file_uniprot = args.output_uniprot
    else:
        output_file_uniprot = None

    if args.output_gencode:
        output_file_gencode = args.output_gencode
    else:
        output_file_gencode = None

    if args.output_canonical:
        output_file_canonical = args.output_canonical
    else:
        output_file_canonical = None
    
    protein_idx_to_seq_dict = get_protein_idx_to_seq_dict(ref_file)
    protein_idx_to_all_tx_dict,protein_orf_type_dict = get_protein_idx_to_all_tx_dict(index_file)
    protein_idx_to_uniprot_dict = get_protein_idx_to_uniprot_dict(index_file)
    protein_seq_dict_uniprot, protein_seq_dict_gencode, protein_seq_dict_canonical, protein_seq_dict_novel = get_separate_protein_seq_dicts(protein_orf_type_dict)
    #print('UniProt proteins in reference: %d'%(len(protein_seq_dict_uniprot)))
    #print('GENCODE proteins in reference: %d'%(len(protein_seq_dict_gencode)))
    #print('Other canonical proteins in reference: %d'%(len(protein_seq_dict_canonical)))
    #print('Novel proteins in reference: %d'%(len(protein_seq_dict_novel)))
    #novel_peptide_count_dict,novel_peptide_mapping_dict = get_checked_peptide_mapping_dicts(checked_file)

    peptide_to_protein_idx_dict,peptide_count_dict,sample_list = get_peptide_mapping_dicts(peptide_file) #, exclude_list=set(novel_peptide_count_dict.keys()))

    if output_file_uniprot:
        n_uniprot1, selected_peptide_dict1, uniprot_peptide_all_mapping_dict1, uniprot_peptide_all_mapping_dict_starts1, uniprot_peptide_all_mapping_dict_ends1 = first_uniprot_filter(peptide_to_protein_idx_dict, protein_orf_type_dict, protein_seq_dict_uniprot, protein_idx_to_uniprot_dict, report_all_hits=True)
        n_uniprot2, selected_peptide_dict2, uniprot_peptide_all_mapping_dict2, uniprot_peptide_all_mapping_dict_starts2, uniprot_peptide_all_mapping_dict_ends2 = second_uniprot_filter(selected_peptide_dict1, protein_seq_dict_uniprot, protein_idx_to_uniprot_dict, report_all_hits=True)
    else:
        n_uniprot1, selected_peptide_dict1 = first_uniprot_filter(peptide_to_protein_idx_dict, protein_orf_type_dict, protein_seq_dict_uniprot, protein_idx_to_uniprot_dict, report_all_hits=False)
        n_uniprot2, selected_peptide_dict2 = second_uniprot_filter(selected_peptide_dict1, protein_seq_dict_uniprot, protein_idx_to_uniprot_dict, report_all_hits=False)
    
    if output_file_gencode:
        n_gencode, selected_peptide_dict3, gencode_peptide_all_mapping_dict, gencode_peptide_all_mapping_dict_starts, gencode_peptide_all_mapping_dict_ends = third_gencode_filter(selected_peptide_dict2, protein_seq_dict_gencode, report_all_hits=True)
    else:
        n_gencode, selected_peptide_dict3 = third_gencode_filter(selected_peptide_dict2, protein_seq_dict_gencode, report_all_hits=False)
    
    if output_file_canonical:
        n_canonical, selected_peptide_dict4, canonical_peptide_all_mapping_dict, canonical_peptide_all_mapping_dict_starts, canonical_peptide_all_mapping_dict_ends = fourth_canonical_filter(selected_peptide_dict3, protein_seq_dict_canonical, report_all_hits=True)
    else:
        n_canonical, selected_peptide_dict4 = fourth_canonical_filter(selected_peptide_dict3, protein_seq_dict_canonical, report_all_hits=False)
 
    n_novel,selected_peptide_dict5, novel_peptide_all_mapping_dict, novel_peptide_all_mapping_dict_starts, novel_peptide_all_mapping_dict_ends = get_all_mappings_of_novel_peptides(selected_peptide_dict4, protein_seq_dict_novel)
    
    print('%s\t%s\t%s\t%s\t%s'%('#Sample', 'UniProt_peptides', 'GENCODE_peptides', 'Other_canonical_peptides', 'novel_peptides'))
    print('%s\t%d\t%d\t%d\t%d'%(peptide_file, n_uniprot1+n_uniprot2, n_gencode, n_canonical, n_novel))
    
    if output_file_novel:
        outf_table(novel_peptide_all_mapping_dict, novel_peptide_all_mapping_dict_starts, novel_peptide_all_mapping_dict_ends, output_file_novel)

    if output_file_canonical:
        outf_table(canonical_peptide_all_mapping_dict, canonical_peptide_all_mapping_dict_starts, canonical_peptide_all_mapping_dict_ends, output_file_canonical)

    if output_file_gencode:
        outf_table(gencode_peptide_all_mapping_dict, gencode_peptide_all_mapping_dict_starts, gencode_peptide_all_mapping_dict_ends, output_file_gencode)

    if output_file_uniprot:
        outf_table(uniprot_peptide_all_mapping_dict1, uniprot_peptide_all_mapping_dict_starts1, uniprot_peptide_all_mapping_dict_ends1, output_file_uniprot)
        outf_table(uniprot_peptide_all_mapping_dict2, uniprot_peptide_all_mapping_dict_starts2, uniprot_peptide_all_mapping_dict_ends2, output_file_uniprot, append=True)
    
