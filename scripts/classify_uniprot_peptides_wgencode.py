import re,os,argparse
import numpy as np

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

def check_shareness(peptide_seq, protein_seq_set):
    if len(protein_seq_set) == 0:
        return True
    for protein_seq in protein_seq_set:
        if protein_seq.find(peptide_seq) < 0:
            return False
    return True

def get_pep_out_line(psm_count_dict, sample_list):
    pep_dict = {}
    for peptide, sample_psm_dict in psm_count_dict.items():
        for sample, psm_count in sample_psm_dict.items():
            if psm_count > 0:
                pep_dict.setdefault(sample, set()).add(peptide)
    pep_list = ['%d'%len(pep_dict.get(sample, set())) for sample in sample_list]
    return '%d\t%s'%(len(psm_count_dict), '\t'.join(pep_list))


def get_psm_out_line(total_psm_dict, sample_list):
    total_psm = sum(total_psm_dict.values())
    psm_list = ['%d'%(total_psm_dict.get(sample, 0)) for sample in sample_list]
    return '%d\t%s'%(total_psm, '\t'.join(psm_list))


def get_length_out_line(total_len_dict, total_psm_dict, sample_list):
    if sum(total_psm_dict.values()) == 0:
        return '0\t%s'%('\t'.join(['0']*len(sample_list)))
    avg_len = sum(total_len_dict.values())/sum(total_psm_dict.values())
    avg_len_list = []
    for sample in sample_list:
        if total_psm_dict.get(sample, 0) == 0:
            avg_len_list.append('0')
        else:
            avg_len_list.append('%2f'%(total_len_dict[sample]/total_psm_dict[sample]))
    return '%.2f\t%s'%(avg_len, '\t'.join(avg_len_list))


def get_len_psm_dicts(psm_count_dict):
    total_len_dict = {}
    total_psm_dict = {}
    for peptide_seq, psm_dict in psm_count_dict.items():
        peptide_len = len(peptide_seq)
        for sample, psm_count in psm_dict.items():
            total_len_dict[sample] = total_len_dict.get(sample, 0) + peptide_len*psm_count
            total_psm_dict[sample] = total_psm_dict.get(sample, 0) + psm_count
    return total_len_dict,total_psm_dict


def outf_table(peptide_list, peptide_mapping_dict, peptide_count_dict, sample_list, output_file=None):
    out_line = '#%s\t%s\t%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id', 'total_count', '\t'.join(sample_list))
    if output_file:
        outf = open(output_file, 'w')
        outf.write(out_line)
    else:
        print(out_line)
    n = 0
    for peptide_seq in peptide_list:
        start, end, tx_id = peptide_mapping_dict[peptide_seq]
        count_dict = peptide_count_dict[peptide_seq]
        count_list = [count_dict[_] for _ in sample_list]
        total_count = sum(count_list)
        count_list = [str(_) for _ in count_list]
        out_line = '%s\t%s\t%s\t%s\t%d\t%s\n'%(peptide_seq, start, end, tx_id, total_count, '\t'.join(count_list))
        if output_file:
            outf.write(out_line)
        else:
            print(out_line)
    if output_file:
        outf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Separate UniProt peptides into three groups.')
    parser.add_argument('-pep', '--peptide_file', help='Input peptide table. Required.', required=True)
    parser.add_argument('-f', '--gencode_protein_fa', help='FASTA of protein sequences in GENCODE. Required', required=True)
    parser.add_argument('-g', '--gencode_gtf', help='GENCODE GTF. Required', required=True)
    parser.add_argument('-u', '--uniprot_fa', help='FASTA of UniProt protein sequences. Required', required=True)
    parser.add_argument('-o1','--output_file1', help='Output peptides mapped to constitutive regions of UniProt and GENOCODE proteins (e.g. shared by all protein isoforms of that gene)', type=str)
    parser.add_argument('-o2','--output_file2', help='Output peptides mapped to alternative regions of UniProt and GENOCODE proteins (e.g. NOT shared by all protein isoforms of that gene)', type=str)
    parser.add_argument('-o3','--output_file3', help='Output peptides mapped to conflict regions of multiple genes.', type=str)
    parser.add_argument('-s','--output_stats', help='Output statistical summary.', type=str)
    args = parser.parse_args()

    peptide_file = args.peptide_file
    gencode_protein_fa = args.gencode_protein_fa
    gencode_gtf = args.gencode_gtf
    uniprot_fa = args.uniprot_fa
    output_file1 = args.output_file1
    output_file2 = args.output_file2
    output_file3 = args.output_file3
    if args.output_stats:
        output_file_stats = args.output_stats
    else:
        output_file_stats = None

    protein_seq_dict = {}
    for line in open(gencode_protein_fa, 'r'):
        if line[0] == '>':
            seq_id = line.rstrip('\n')[1:]
        else:
            protein_seq_dict[seq_id] = protein_seq_dict.get(seq_id, '') + line.rstrip('\n').upper()

    tx_to_gene_name_dict = {}
    for line in open(gencode_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                #print('Excluded line in annotation GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            gene_name_list = d.get('gene_name', [''])
            tx_to_gene_name_dict[tx_id] = gene_name_list

    gene_protein_seq_dict = {}
    for tx_id, protein_seq in protein_seq_dict.items():
        gene_name_list = tx_to_gene_name_dict[tx_id]
        for gene_name in gene_name_list:
            gene_protein_seq_dict.setdefault(gene_name, set()).add(protein_seq)

    uniprot_protein_seq_dict = {}
    uniprot_gene_name_dict = {}
    for line in open(uniprot_fa, 'r'):
        if line[0] == '>':
            seq_id = line.rstrip('\n').split('|')[1]
            match = re.findall('GN=(\S*)', line)
            if len(match) == 0:
                gene_name = ''
            else:
                gene_name = match[0]
            uniprot_gene_name_dict[seq_id] = gene_name
        else:
            if gene_name != '':
                uniprot_protein_seq_dict[seq_id] = uniprot_protein_seq_dict.get(seq_id, '') + line.rstrip('\n')

    peptide_gene_name_dict = {}
    peptide_mapping_dict = {}
    peptide_count_dict = {}
    sample_list = []
    out_line = '#%s\t%s\t%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id', 'total_count', '\t'.join(sample_list))
    header_dict = {}
    for line in open(peptide_file, 'r'):
        if line.startswith('#'):
            arr = line.rstrip('\n')[1:].split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
                if arr[i] not in ['peptide', 'starts', 'ends', 'tx_id', 'total_count']:
                    sample_list.append(arr[i])
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[0]
        peptide_mapping_dict[peptide_seq] = (arr[header_dict['starts']], arr[header_dict['ends']], arr[header_dict['tx_id']])
        seq_id_list = arr[header_dict['tx_id']].split(',')
        for seq_id in seq_id_list:
            gene_name = uniprot_gene_name_dict[seq_id]
            peptide_gene_name_dict.setdefault(peptide_seq, set()).add(gene_name)
        for sample in sample_list:
            peptide_count_dict.setdefault(peptide_seq, {})[sample] = int(arr[header_dict[sample]])

    shared_peptides = set()
    nonshared_peptides = set()
    conflict_peptides = set()
    for peptide_seq, gene_set in peptide_gene_name_dict.items():
        shareness_res = set()
        for gene_name in gene_set:
            gencode_protein_seq_set = gene_protein_seq_dict.get(gene_name, set())
            shareness_res.add(check_shareness(peptide_seq, gencode_protein_seq_set))
        if len(shareness_res) > 1:
            conflict_peptides.add(peptide_seq)
        elif shareness_res == set([True]):
            shared_peptides.add(peptide_seq)
        else:
            nonshared_peptides.add(peptide_seq)

    outf_table(shared_peptides, peptide_mapping_dict, peptide_count_dict, sample_list, output_file1)
    outf_table(nonshared_peptides, peptide_mapping_dict, peptide_count_dict, sample_list, output_file2)
    outf_table(conflict_peptides, peptide_mapping_dict, peptide_count_dict, sample_list, output_file3)

    shared_peptide_count_dict = {k:peptide_count_dict[k] for k in shared_peptides}
    nonshared_peptide_count_dict = {k:peptide_count_dict[k] for k in nonshared_peptides}
    conflict_peptide_count_dict = {k:peptide_count_dict[k] for k in conflict_peptides}

    total_len_dict_shared, total_psm_dict_shared = get_len_psm_dicts(shared_peptide_count_dict)
    total_len_dict_nonshared, total_psm_dict_nonshared = get_len_psm_dicts(nonshared_peptide_count_dict)
    total_len_dict_conflict, total_psm_dict_conflict = get_len_psm_dicts(conflict_peptide_count_dict)

    if output_file_stats:
        output = open(output_file_stats, 'w')
        output.write('## Number of samples: %d\n'%(len(sample_list)))
        output.write('\n')
        output.write('## Number of peptides in all/each sample:\n')
        output.write('#Peptides_types\tPeptides_count\t%s\n'%('\t'.join(sample_list)))
        out_line = get_pep_out_line(shared_peptide_count_dict, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Shared)', out_line))

        out_line = get_pep_out_line(nonshared_peptide_count_dict, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Nonshared)', out_line))

        out_line = get_pep_out_line(conflict_peptide_count_dict, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Conflict_Genes)', out_line))

        output.write('\n')
        output.write('## Number of PSMs in all/each sample:\n')
        output.write('#Peptides_types\tTotal_PSM\t%s\n'%('\t'.join(sample_list)))
        out_line = get_psm_out_line(total_psm_dict_shared, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Shared)', out_line))

        out_line = get_psm_out_line(total_psm_dict_nonshared, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Nonshared)', out_line))

        out_line = get_psm_out_line(total_psm_dict_conflict, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Conflict_Genes)', out_line))
        output.write('\n')
        output.write('## Average PSM length in all/each sample:\n')
        output.write('#Peptides_types\tAvg_len\t%s\n'%('\t'.join(sample_list)))
        out_line = get_length_out_line(total_len_dict_shared, total_psm_dict_shared, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Shared)', out_line))

        out_line = get_length_out_line(total_len_dict_nonshared, total_psm_dict_nonshared, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Nonshared)', out_line))

        out_line = get_length_out_line(total_len_dict_conflict, total_psm_dict_conflict, sample_list)
        output.write('%s\t%s\n'%('UniProt&GENCODE(Conflict_Genes)', out_line))
        output.close()

        


