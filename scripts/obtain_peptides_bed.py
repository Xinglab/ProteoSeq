import os,sys,argparse
import re

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

def genome_position(region_start, region_end, strand, exon_list):
    if region_start >= region_end:
        print('Wrong input region coordinates!')
        return (0, 0, [])
    exon_list = sorted(exon_list)
    genome_start, genome_end = region_start, region_end
    split_region_list = []
    if len(exon_list) <= 1:
        exon_start, exon_end = exon_list[0]
        if (exon_end-exon_start+1) < region_end:
            print('Region exceed the total length of exons!')
            return (0, 0, [])
        if strand == '+':
            genome_start = exon_start + region_start - 1
            genome_end = exon_start + region_end - 1
            split_region_list = [(genome_start, genome_end)]
        elif strand == '-':
            genome_start = exon_end - region_start + 1
            genome_end = exon_end - region_end + 1
            split_region_list = [(genome_end, genome_start)]
        return (genome_start, genome_end, split_region_list)
    checked_len = 0 
    start_record_tag = 0
    end_record_tag = 0
    if strand == '+':
        for (exon_start, exon_end) in exon_list:
            exon_len = exon_end - exon_start + 1
            if checked_len < region_start <= checked_len+exon_len:  #in this region
                genome_start = exon_start + region_start - checked_len - 1
                start_record_tag = 1
            if checked_len < region_end <= checked_len+exon_len:
                genome_end = exon_start + region_end - checked_len - 1
                end_record_tag = 1
            checked_len = checked_len + exon_len
            # split_record_tag: 1: find start position; 2: find stop position; 0.5: in the middle part
            if (start_record_tag==1) and (end_record_tag==0):  
                split_region = (genome_start, exon_end)
                start_record_tag = 2
            elif (start_record_tag==1) and (end_record_tag == 1):
                split_region = (genome_start, genome_end)
                start_record_tag = 2
                end_record_tag = 2
            elif (start_record_tag==2) and (end_record_tag == 0):
                split_region = (exon_start, exon_end)
            elif (start_record_tag==2) and (end_record_tag == 1):
                split_region = (exon_start, genome_end)
                end_record_tag = 2
            if start_record_tag == 2:
                split_region_list.append(split_region)
                if end_record_tag==2:
                    break
    elif strand == '-':
        first_start = exon_list[0][0]
        second_start = exon_list[1][0]
        if first_start < second_start:
            exon_list_neg = exon_list[::-1]
        else:
            exon_list_neg = exon_list
        for (exon_start, exon_end) in exon_list_neg:
            exon_len = exon_end - exon_start + 1
            if checked_len < region_start <= checked_len+exon_len:  #in this region
                genome_start = exon_end - region_start + checked_len + 1
                start_record_tag = 1
            if checked_len < region_end <= checked_len+exon_len:
                genome_end = exon_end - region_end + checked_len + 1 #compatible with gtf coordinate
                end_record_tag = 1
            checked_len = checked_len + exon_len
            # split_record_tag: 1: find start position; 2: find stop position; 0.5: in the middle part
            if (start_record_tag==1) and (end_record_tag==0):
                split_region = (exon_start, genome_start)
                start_record_tag = 2
            elif (start_record_tag==1) and (end_record_tag == 1):
                split_region = (genome_end, genome_start)
                start_record_tag = 2
                end_record_tag = 2
            elif (start_record_tag==2) and (end_record_tag == 0):
                split_region = (exon_start, exon_end)
            elif (start_record_tag==2) and (end_record_tag == 1):
                split_region = (genome_end, exon_end)
                end_record_tag = 2
            if start_record_tag == 2:
                split_region_list.append(split_region)
                if end_record_tag==2:
                    split_region_list = split_region_list[::-1] #from smaller coords to larger
                    break
    if (start_record_tag != 2) or (end_record_tag != 2):
        print('Region exceed the total length of exons!')
        return (0, 0, [])
    return (genome_start, genome_end, split_region_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate index table with information from GTF.')
    parser.add_argument('-i', '--input_table', help='Input table of peptides. Required', type=str, required=True)
    parser.add_argument('-f', '--protein_fasta', help='Input FASTA of protein sequences. Required', type=str, required=True)
    parser.add_argument('-idx', '--index_file', help='Input file mapping protein index to of protein sequence. Required', type=str, required=True)
    parser.add_argument('-g', '--annotation_gtf', help='Input GTF for annotation. Required', type=str, required=True)
    parser.add_argument('-o', '--output_bed', help='Output BED file of transcripts with peptides evidence. Required', type=str, required=True)
    args = parser.parse_args()
    
    input_table = args.input_table
    annotation_gtf = args.annotation_gtf
    protein_fasta = args.protein_fasta
    index_file = args.index_file
    output_bed = args.output_bed

    protein_idx_seq_dict = {}
    for line in open(protein_fasta, 'r'):
        if line[0]=='>':
            protein_idx = int(line.rstrip('\n').replace('>Protein_', ''))
        else:
            protein_idx_seq_dict[protein_idx] = protein_idx_seq_dict.get(protein_idx, '') + line.rstrip('\n')

    protein_idx_all_tx_dict = {}
    header_dict = {}
    for line in open(index_file, 'r'):
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = int(arr[header_dict['protein_idx']].replace('Protein_', ''))
        tx_id = arr[header_dict['transcript_id']]
        if tx_id == '':
            continue
        protein_idx_all_tx_dict.setdefault(protein_idx, set()).add(tx_id)

    tx_coord_dict = {}
    cds_coord_dict = {}
    for line in open(annotation_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            tx_coord_dict[tx_id] = (arr[0], arr[6], int(arr[3]), int(arr[4]))
        if arr[2] == 'CDS':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            cds_coord_dict.setdefault(tx_id, []).append((int(arr[3]), int(arr[4])))

    protein_all_mappings_dict = {}
    header_dict = {}
    for line in open(input_table, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[0]
        protein_idx_list = arr[1].split(';')
        for protein_idx in protein_idx_list:
            protein_idx = int(protein_idx)
            protein_seq = protein_idx_seq_dict[protein_idx]
            for match in re.finditer(peptide_seq, protein_seq):
                start = match.start()+1 #Python mode
                end = match.end()
                protein_all_mappings_dict.setdefault(protein_idx, {})[(start, end)] = peptide_seq

    n = 0
    tx_all_mappings_dict = {}
    no_expressed_tx = set()
    for protein_idx in sorted(protein_all_mappings_dict):
        temp_d = protein_all_mappings_dict[protein_idx]
        tx_id_list = protein_idx_all_tx_dict.get(protein_idx, [])
        if len(tx_id_list) == 0:
            no_expressed_tx.add(protein_idx)
            continue
        for tx_id in tx_id_list:
            if len(tx_all_mappings_dict.get(tx_id, {})):
                print(tx_id)
                break
            tx_all_mappings_dict[tx_id] = temp_d
        n += 1

    print('No expressed tx:%d, Expressed tx:%d'%(len(no_expressed_tx), n))

    peptide_bed_all_idx_dict = {}
    for tx_id, temp_d in tx_all_mappings_dict.items():
        for peptide_coord in sorted(temp_d):
            peptide_start, peptide_end = peptide_coord
            peptide_seq = temp_d[peptide_coord]
            peptide_start_in_cds = (peptide_start-1)*3 + 1
            peptide_end_in_cds = (peptide_end)*3
            (chrom, strand, tx_start, tx_end) = tx_coord_dict[tx_id]
            cds_coord_list = cds_coord_dict[tx_id]
            peptide_start_genome, peptide_end_genome, peptide_exon_coords = genome_position(peptide_start_in_cds, peptide_end_in_cds, strand, cds_coord_list)
            if peptide_start_genome > peptide_end_genome:
                peptide_start_genome, peptide_end_genome = peptide_end_genome, peptide_start_genome
            peptide_genomic_coord = (chrom, peptide_start_genome, peptide_end_genome, strand)
            peptide_bed_all_idx_dict.setdefault(peptide_genomic_coord, {}).setdefault(peptide_seq, set()).add(tx_id)

    output = open(output_bed, 'w')
    for peptide_genomic_coord in sorted(peptide_bed_all_idx_dict):
        (chrom, peptide_start_genome, peptide_end_genome, strand) = peptide_genomic_coord
        peptide_all_tx_dict = peptide_bed_all_idx_dict[peptide_genomic_coord]
        for peptide_seq, tx_list in peptide_all_tx_dict.items():
            n = output.write('%s\t%d\t%d\t%s_%s\t0\t%s\n'%(chrom, peptide_start_genome-1, peptide_end_genome, peptide_seq, ';'.join(tx_list), strand))
    output.close()



