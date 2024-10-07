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


def get_separate_protein_seq_dicts(protein_idx_to_seq_dict, protein_idx_to_uniprot_dict):
    protein_seq_dict_canonical = {}
    protein_seq_dict_isoform = {}
    for protein_idx, tx_id_list in protein_idx_to_uniprot_dict.items():
        protein_seq = protein_idx_to_seq_dict[protein_idx]
        for tx_id in tx_id_list:
            if tx_id.find('-') >= 0:
                protein_seq_dict_isoform[tx_id] = protein_seq
            else:
                protein_seq_dict_canonical[tx_id] = protein_seq
    return protein_seq_dict_canonical, protein_seq_dict_isoform


def get_peptide_mapping_dicts(input_file):
    peptide_mapping_dict = {}
    peptide_abun_dict = {}
    sample_list = []
    header_dict = {}
    header_rev_dict = {}
    abun_idx = 1
    for line in open(input_file, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
                header_rev_dict[i] = arr[i]
                if arr[i].startswith('total_'):
                    if abun_idx > 1:
                        continue
                    abun_idx = i
            sample_list = arr[abun_idx:]
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[header_dict['peptide']]
        peptide_abun_dict[peptide_seq] = arr[abun_idx:]
        peptide_mapping_dict[peptide_seq] = (arr[header_dict['starts']], arr[header_dict['ends']], arr[header_dict['tx_id']])
    return peptide_mapping_dict,peptide_abun_dict,sample_list


def outf_table(peptide_list, peptide_mapping_dict, peptide_abun_dict, sample_list, output_file):
    out_line = '#%s\t%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id', '\t'.join(sample_list))
    outf = open(output_file, 'w')
    outf.write(out_line)
    for peptide_seq in peptide_list:
        start, end, tx_id = peptide_mapping_dict[peptide_seq]
        count_list = peptide_abun_dict[peptide_seq]
        out_line = '%s\t%s\t%s\t%s\t%s\n'%(peptide_seq, start, end, tx_id, '\t'.join(count_list))
        outf.write(out_line)
    outf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Separate UniProt peptides into three groups.')
    parser.add_argument('-pep', '--peptide_file', help='Input peptide table. Required.', required=True)
    parser.add_argument('-ref', '--ref_file', help='Input FASTA file of UniProt reference protein sequences. Required', required=True)
    parser.add_argument('-idx', '--index_file', help='Input file mapping indexes to reference protein sequences . Required', required=True)
    parser.add_argument('-o1','--output_file1', help='Output peptides mapped to canonical UniProt proteins', type=str)
    parser.add_argument('-o2','--output_file2', help='Output peptides uniquely mapped to UniProt alternative protein isoforms', type=str)
    args = parser.parse_args()

    index_file = args.index_file
    ref_file = args.ref_file
    peptide_file = args.peptide_file
    output_file1 = args.output_file1
    output_file2 = args.output_file2

    protein_idx_to_seq_dict = get_protein_idx_to_seq_dict(ref_file)
    protein_idx_to_uniprot_dict = get_protein_idx_to_uniprot_dict(index_file)
    protein_seq_dict_canonical,protein_seq_dict_isoform = get_separate_protein_seq_dicts(protein_idx_to_seq_dict, protein_idx_to_uniprot_dict)

    protein_canonical_id_to_isoform_seq_dict = {}
    for seq_id, seq in protein_seq_dict_isoform.items():
        canonical_id = seq_id.split('-')[0]
        protein_canonical_id_to_isoform_seq_dict.setdefault(canonical_id, set()).add(seq)

    peptide_mapping_dict, peptide_abun_dict, sample_list = get_peptide_mapping_dicts(peptide_file)

    canonical_peptides = set()
    isoform_peptides = set()
    for peptide_seq, peptide_mapping in peptide_mapping_dict.items():
        for protein_id, protein_seq in protein_seq_dict_canonical.items():
            start = protein_seq.find(peptide_seq) + 1
            if start >= 1:
                canonical_peptides.add(peptide_seq)
                break
        else:
            isoform_peptides.add(peptide_seq)

    outf_table(canonical_peptides, peptide_mapping_dict, peptide_abun_dict, sample_list, output_file1)
    outf_table(isoform_peptides, peptide_mapping_dict, peptide_abun_dict, sample_list, output_file2)

