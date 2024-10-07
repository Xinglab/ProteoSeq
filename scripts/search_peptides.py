import os, argparse, re
from Bio.Seq import Seq
from Bio import SeqIO

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

def fetch_seq(genome, chrom, start, end, strand):
    seq = genome[chrom][start:end]
    if strand == "-":
        seq = seq.reverse_complement()
    return str(seq.seq).upper()

def get_seq_from_exon_list(chrom, strand, exon_list):
    exon_list = sorted(exon_list)
    final_seq = ''
    if len(exon_list) <= 1:
        exon_start, exon_end = exon_list[0]
        if exon_start > exon_end:
            return '',[]
        final_seq = fetch_seq(genome, chrom, exon_start-1, exon_end, strand)
        return final_seq
    if strand == '-':
        exon_list = exon_list[::-1] #reversed
    for (exon_start,exon_end) in exon_list:
        if exon_start > exon_end:
            return '',[]
        temp_seq = fetch_seq(genome, chrom, exon_start-1, exon_end, strand)
        final_seq += temp_seq
    return final_seq

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

def check_overlap(coord1, coord2):
    sorted_coord = sorted([coord1[0], coord1[1], coord2[0], coord2[1]])
    if sorted_coord == [coord1[0], coord1[1], coord2[0], coord2[1]]:
        if coord1[1] == coord2[0]: #One nucleotide overlapping
            return True
        return False
    if sorted_coord == [coord2[0], coord2[1], coord1[0], coord1[1]]:
        if coord2[1] == coord1[0]: #One nucleotide overlapping
            return True
        return False
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search peptides in a given GTF files.')
    parser.add_argument('-i', '--input_fasta', help='Input fasta of peptide sequences. Required', type=str, required=True)
    parser.add_argument('-c', '--cds_gtf', help='Input GTF of searching space with CDS coordinates', type=str, required=True)
    parser.add_argument('-g', '--genome_fasta', help='Input fasta of the genome sequence. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF with peptide mappings. Required', type=str, required=True)
    parser.add_argument('--output_table', help='Output tabular table with peptide mappings', type=str)
    parser.add_argument('--output_unmapped', help='Output FA of peptides not mapped.')
    parser.add_argument('-v', '--verbose', help='Set to zero to silent stdout skipped lines, default=0', type=int, choices=[0, 1], default=0)
    args = parser.parse_args()
    
    input_fasta = args.input_fasta
    cds_gtf = args.cds_gtf
    genome_fasta = args.genome_fasta
    output_file_gtf = args.output_gtf
    verbose = args.verbose
    if args.output_table:
        output_file_table = args.output_table
    else:
        output_file_table = None
    if args.output_unmapped:
        output_file_unmapped = args.output_unmapped
    else:
        output_file_unmapped = None
    
    peptide_seq_dict = {}
    for line in open(input_fasta, 'r'):
        if line.startswith('>'):
            peptide_id = line.rstrip('\n')[1:]
        else:
            peptide_seq_dict[peptide_id] = peptide_seq_dict.get(peptide_id, '') + line.rstrip('\n')

    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    tx_cds_coord_dict = {}
    tx_exon_coord_to_anno_dict = {}
    tx_chrom_dict = {}
    tx_strand_dict = {}
    tx_anno_dict = {}
    for line in open(cds_gtf,'r'):
        if line[0] == '#':
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('No transcript ID in line in gtf: %s'%line)
                continue
            tx_anno_dict[tx_id] = arr[8]
        elif arr[2] == 'CDS':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('No transcript ID in line in gtf: %s'%line)
                continue
            cds_coord = (int(arr[3]), int(arr[4]))
            tx_cds_coord_dict.setdefault(tx_id, []).append(cds_coord)
            tx_chrom_dict[tx_id] = arr[0]
            tx_strand_dict[tx_id] = arr[6]
        elif arr[2] == 'exon':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('No transcript ID in line in gtf: %s'%line)
                continue
            exon_coord = (int(arr[3]), int(arr[4]))
            tx_exon_coord_to_anno_dict.setdefault(tx_id, {})[exon_coord] = arr[8]

    tx_cds_coord_dict = {k:sorted(v) for k,v in tx_cds_coord_dict.items()}
    tx_list_wCDS = sorted(tx_cds_coord_dict)
    if len(tx_list_wCDS) == 0:
        exit('No CDS is detected in the GTF file! Make sure the 2nd column of GTF has CDS')

    protein_seq_dict = {}
    for tx_id in tx_list_wCDS:
        chrom = tx_chrom_dict[tx_id]
        strand = tx_strand_dict[tx_id]
        cds_coord_list = tx_cds_coord_dict.get(tx_id, [])
        if len(cds_coord_list) == 0: #No CDS is found in annotation
            continue
        cds_seq = get_seq_from_exon_list(chrom, strand, cds_coord_list)
        if len(cds_seq) % 3 > 0:
            if verbose == 1:
                print('%s does not have a valid CDS sequence, length is not multiple of three, ignored.'%tx_id)
            continue
        protein_seq=Seq(cds_seq).translate()
        if (protein_seq.find('*') >= 0) and (protein_seq.find('*') < (len(protein_seq)-1)):
            if verbose == 1:
                print('%s CDS sequence has in-frame stop codons, ignored.'%tx_id)
            continue
        protein_seq_dict[tx_id] = str(protein_seq).upper()

    protein_all_tx_dict = {}
    for tx_id, protein_seq in protein_seq_dict.items():
        protein_all_tx_dict.setdefault(protein_seq, []).append(tx_id)

    peptide_mapping_dict = {}
    for peptide_id, peptide_seq in peptide_seq_dict.items():
        for protein_seq, all_tx_list in protein_all_tx_dict.items():
            for match in re.finditer(peptide_seq, protein_seq):
                start = match.start()+1 #Python mode
                end = match.end()
                for tx_id in all_tx_list: 
                    peptide_mapping_dict.setdefault(peptide_id, {}).setdefault(tx_id, []).append((start, end))

    if output_file_unmapped:
        output_unmapped = open(output_file_unmapped, 'w')
        unmapped_peptide = set(peptide_seq_dict) - set(peptide_mapping_dict)
        for peptide_id in sorted(unmapped_peptide):
            peptide_seq = peptide_seq_dict[peptide_id]
            output_unmapped.write('>%s\n%s\n'%(peptide_id, peptide_seq))
        output_unmapped.close()

    if output_file_table:
        output_table = open(output_file_table, 'w')
        output_table.write('#%s\t%s\t%s\t%s\n'%('peptide', 'starts', 'ends', 'tx_id'))
    else:
        output_table = None

    output = open(output_file_gtf, 'w')
    for peptide_id, mapping_dict in peptide_mapping_dict.items():
        peptide_seq = peptide_seq_dict[peptide_id]
        all_tx_list = []
        all_starts_list = []
        all_ends_list = []
        for tx_id, mapping_coord_list in mapping_dict.items():
            chrom = tx_chrom_dict[tx_id]
            strand = tx_strand_dict[tx_id]
            attribute_output = ['peptide_id "%s"'%peptide_id, 'peptide_seq "%s"'%peptide_seq]
            tx_anno = tx_anno_dict[tx_id]
            attribute_output.append(tx_anno)
            attribute = '; '.join(attribute_output)
            cds_coord_list = tx_cds_coord_dict[tx_id]
            for mapping_coord in mapping_coord_list:
                all_tx_list.append(tx_id)
                all_starts_list.append(str(mapping_coord[0]))
                all_ends_list.append(str(mapping_coord[1]))
                peptide_start_in_cds = (mapping_coord[0]-1)*3 + 1
                peptide_end_in_cds = (mapping_coord[1])*3
                peptide_start_genome, peptide_end_genome, peptide_exon_coords = genome_position(peptide_start_in_cds, peptide_end_in_cds, strand, cds_coord_list)
                for peptide_coord in peptide_exon_coords:
                    exon_coord_anno_dict = tx_exon_coord_to_anno_dict.get(tx_id, {})
                    for exon_coord in sorted(exon_coord_anno_dict):
                        if check_overlap(peptide_coord, exon_coord):
                            exon_anno = exon_coord_anno_dict[exon_coord]
                            attribute_output_exon = ['peptide_id "%s"'%peptide_id, 'peptide_seq "%s"'%peptide_seq]
                            attribute_output_exon.append(exon_anno)
                            attribute_exon = '; '.join(attribute_output_exon)
                            gtf_peptide = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n'%(chrom, 'ProteoSeq', 'peptide', peptide_coord[0], peptide_coord[1], '.', strand, '.', attribute_exon)
                            output.write(gtf_peptide)
                            break
                    else:
                        gtf_peptide = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n'%(chrom, 'ProteoSeq', 'peptide', peptide_coord[0], peptide_coord[1], '.', strand, '.', attribute)
                        output.write(gtf_peptide)
                        if verbose == 1:
                            print('Peptide-%s %s in %s (%d, %d) does not have any overlaping CDS exons'%(peptide_id, peptide_seq, tx_id, peptide_coord[0], peptide_coord[1]))
        if output_table:
            output_table.write('%s\t%s\t%s\t%s\n'%(peptide_seq, ','.join(all_starts_list), ','.join(all_ends_list), ','.join(all_tx_list)))
    output.close()

    if output_table:
        output_table.close()

