import os,sys,argparse
import re
from Bio.Seq import Seq
from Bio import SeqIO

def parse_attributes(string):
    d = {}
    for i in string.split(';'):
        if i == '':
            continue
        i_list = i.strip().split(' ')
        if len(i_list) < 2:
            #print(string)
            continue
        d.setdefault(i_list[0], []).append(i_list[1].strip('"'))
    return d

def filter_attributes(d, inc_list=[]):
    new_d = {}
    for k, v in d.items():
        if k in inc_list:
            new_d[k] = v
    return new_d

def deparse_attributes(d):
    string = []
    for k,v_list in d.items():
        for v in v_list:
            string.append('%s \"%s\"'%(k, v))
    return '; '.join(string)+';'

def fetch_seq(genome, chrom, start, end, strand):
    seq = genome[chrom][start:end]
    if strand == "-":
        seq = seq.reverse_complement()
    return str(seq.seq).upper()

def get_seq_from_exon_list(chrom, strand, exon_list, add_tail=0):
    if len(exon_list) == 0:
        return '', []
    exon_list = sorted(exon_list)
    final_seq = ''
    exon_list_wtail = []
    if len(exon_list) <= 1:
        exon_start, exon_end = exon_list[0]
        if exon_start > exon_end:
            return '', []
        final_seq = fetch_seq(genome, chrom, exon_start-1, exon_end, strand)
        if add_tail > 0:
            if strand == '+':
                tail_start = exon_end + 1
                tail = fetch_seq(genome, chrom, tail_start-1, tail_start-1+add_tail, strand)
                exon_list_wtail = [(exon_start, exon_end+add_tail)]
            elif strand == '-':
                tail_start = exon_start
                tail = fetch_seq(genome, chrom, tail_start-1-add_tail, tail_start-1, strand)
                exon_list_wtail = [(exon_start-add_tail, exon_end)]
            final_seq += tail
        return final_seq, exon_list_wtail
    if strand == '-':
        exon_list = exon_list[::-1] #reversed
    exon_list_wtail = exon_list[:] #reversed
    for (exon_start,exon_end) in exon_list:
        if exon_start > exon_end:
            return '', []
        temp_seq = fetch_seq(genome, chrom, exon_start-1, exon_end, strand)
        seq = temp_seq[:-1] + temp_seq[-1].lower()
        final_seq += seq
    #if add_tail <= 0:
    temp_final_seq = final_seq[:-1] + final_seq[-1].upper()
    final_seq = temp_final_seq
    if add_tail <= 0:
        return final_seq, exon_list_wtail
    if strand == '+':
        tail_start = exon_end+1
        tail = fetch_seq(genome, chrom, tail_start-1, tail_start-1+add_tail, strand)
        exon_list_wtail[-1] = (exon_start, exon_end+add_tail)
    elif strand == '-':
        tail_start = exon_start
        tail = fetch_seq(genome, chrom, tail_start-1-add_tail, tail_start-1, strand)
        exon_list_wtail[-1] = (exon_start-add_tail, exon_end) #reversed
        exon_list_wtail = exon_list_wtail[::-1] #reverse back to smaller coord to larger coord
    final_seq += tail
    return final_seq, exon_list_wtail


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

def lower_ss(seq, pos_s0_list):
    if len(pos_s0_list)==0:
        return seq
    for pos_s0 in pos_s0_list:
        if 0 <= pos_s0 < len(seq):
            temp_seq = seq[:pos_s0] + seq[pos_s0].lower() + seq[pos_s0+1:]
            seq = temp_seq
    return seq

def check_NMD_by_seq(tx_seq_w3nt, cds_seq):
    ss_pos_s0 = [i.start() for i in re.finditer("[a-z]", tx_seq_w3nt)]
    ss_pos_s0_in_cds = [i.start() for i in re.finditer("[a-z]", cds_seq)]
    if len(ss_pos_s0) == 0: #Only one exon
        return False
    if len(cds_seq) < 150: # or 200, too short transcript doesn't have NMD
        return False
    last_EJC = ss_pos_s0[-1] + 1
    start_codon = tx_seq_w3nt.upper().find(cds_seq.upper()) + 1
    stop_codon = start_codon + len(cds_seq) - 1
    dist = last_EJC - stop_codon
    if dist < 50: #No NMD
        return False
    return True

def get_longest_ORF_woPTC(start_s0_list, stop_s0_list):
    dist_pairs_dict = {}
    for s0 in sorted(start_s0_list):
        for e0 in sorted(stop_s0_list):
            dist = e0 - s0
            if dist < 0:
                continue
            if dist % 3 > 0:
                continue
            if dist < 60:
                break # shorter than 20AA, discarded
            dist_pairs_dict.setdefault(dist, set()).add((s0, e0))
            break #larger stop would involve the previous PTC
    if len(dist_pairs_dict) == 0:
        return []
    max_dist = max(dist_pairs_dict.keys())
    pair_codon_list = list(sorted(dist_pairs_dict[max_dist]))
    return pair_codon_list

def translate_tx(tx_seq_w3nt, strand, exon_list):
    start_codon_pos_s0_list = [i.start() for i in re.finditer('ATG', tx_seq_w3nt, flags=re.IGNORECASE)]
    stop_codon_pos_s0_list = [i.start() for i in re.finditer("TAG|TGA|TAA", tx_seq_w3nt, flags=re.IGNORECASE)]
    pair_codon_list = get_longest_ORF_woPTC(start_codon_pos_s0_list, stop_codon_pos_s0_list)
    longest_protein_seq = ''
    longest_cds_seq = ''
    longest_cds_exon_coord_list = []
    longest_cds_coord_in_tx = (0, 0)
    if len(pair_codon_list) == 0:
        return longest_protein_seq, longest_cds_seq, longest_cds_exon_coord_list, longest_cds_coord_in_tx
    longest_len = 0
    for pair_codon in pair_codon_list: #loop to get longest orf
        s0, e0 = pair_codon[0], pair_codon[1]
        cds_seq = tx_seq_w3nt[s0:e0+3]
        ss_pos_s0_in_cds = [i.start() for i in re.finditer("[a-z]", cds_seq)]
        cds_start_genome, cds_end_genome, cds_exon_coord_list = genome_position(s0+1, e0+3, strand, exon_list)
        protein_seq = str(Seq(cds_seq).translate())
        ss_pos_s0_in_aa = [_//3 for _ in ss_pos_s0_in_cds]
        protein_seq_with_ss = lower_ss(protein_seq, ss_pos_s0_in_aa)
        if protein_seq_with_ss[-1]=='*': #remove final stop codon-derived AA
            protein_seq_with_ss = protein_seq_with_ss[0:-1]
        protein_len = len(protein_seq_with_ss)
        if protein_len >= longest_len:
            if protein_len == longest_len:
                # if two orf have same length, keep the one with more exons
                if len(cds_exon_coord_list) <= len(longest_cds_exon_coord_list):
                    continue
            longest_protein_seq = protein_seq_with_ss
            longest_cds_seq = cds_seq
            longest_cds_exon_coord_list = cds_exon_coord_list
            longest_cds_coord_in_tx = (s0+1, e0+3)
            longest_len = protein_len
    return longest_protein_seq, longest_cds_seq, longest_cds_exon_coord_list, longest_cds_coord_in_tx


def compare_exon_coords_vs_cds_coords(exon_coord_list, cds_exon_coord_list, strand):
    exon_coord_list = sorted(exon_coord_list)
    cds_exon_coord_list = sorted(cds_exon_coord_list)
    exon_coord_list_w3nt = []
    if strand == '+':
        exon_coord_list_w3nt = exon_coord_list[:-1]
        last_exon = exon_coord_list[-1]
        last_cds_exon = cds_exon_coord_list[-1]
        if last_cds_exon[1] > last_exon[1]:
            exon_coord_list_w3nt.append((last_exon[0], last_cds_exon[1]))
        else:
            exon_coord_list_w3nt.append((last_exon[0], last_exon[1]))
        return exon_coord_list_w3nt
    else:
        exon_coord_list_w3nt = exon_coord_list[1:][::-1]
        last_exon = exon_coord_list[0]
        last_cds_exon = cds_exon_coord_list[0]
        if last_cds_exon[0] < last_exon[0]:
            exon_coord_list_w3nt.append((last_cds_exon[0], last_exon[1]))
        else:
            exon_coord_list_w3nt.append((last_exon[0], last_exon[1]))
        return exon_coord_list_w3nt[::-1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict longest ORF in the transcripts and translate to protein sequences.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file. Attributes in the GTF will be kept in the output GTF. Required', type=str)
    parser.add_argument('-g', '--genome_fasta', help='Input genome fasta file. Required', type=str, required=True)
    parser.add_argument('-o', '--output_gtf', help='Output GTF file of CDS coordinates. Required', type=str, required=True)
    parser.add_argument('--output_protein', help='Output additional fasta file of protein sequences', type=str)
    parser.add_argument('--output_transcript', help='Output additional fasta file of transcripts sequences', type=str)
    parser.add_argument('--output_cds', help='Output additional fasta file of CDS sequences', type=str)
    parser.add_argument('--add_tail', help='Add X-nt tail to the end of transcripts when predict ORF, default=3', type=int, default=3)
    parser.add_argument('--keep_nmd', help='Set to True to output sequences for NMD transcripts. Set to \"Add\" to include NMD sequences (marked in sequence ID) to protein coding sequences.', default='False', choices=['Add', 'True', 'False'])
    parser.add_argument('-v', '--verbose', help='Set to zero to silent stdout skipped lines, default=0', type=int, choices=[0, 1], default=0)
    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    genome_fasta = args.genome_fasta
    add_tail = int(args.add_tail)

    output_gtf = args.output_gtf

    if args.output_protein:
        output_protein = args.output_protein
    else:
        output_protein = None

    if args.output_transcript:
        output_transcript = args.output_transcript
    else:
        output_transcript = None

    if args.output_cds:
        output_cds = args.output_cds
    else:
        output_cds = None

    keep_nmd = args.keep_nmd 
    
    verbose = args.verbose
    
    attribute_list = ['transcript_id', 'transcript_name', 'transcript_type', 'gene_id', 'gene_type', 'gene_name', 'tag']    
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    tx_exon_coord_dict = {}
    chrom_dict = {}
    strand_dict = {}
    tx_attribute_dict = {}
    tx_gtf_dict = {}
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'exon':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('Excluded line in input GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            tx_attribute_dict[tx_id] = d
            tx_exon_coord_dict.setdefault(tx_id, []).append((int(arr[3]), int(arr[4])))
            chrom_dict[tx_id] = arr[0]
            strand_dict[tx_id] = arr[6]
        elif arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('Excluded line in input GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            tx_gtf_dict[tx_id] = line

    tx_exon_coord_dict = {k:sorted(v) for k,v in tx_exon_coord_dict.items()}

    tx_list = sorted(tx_exon_coord_dict)

    tx_list_wCDS = []
    tx_seq_w3nt_dict = {}
    tx_exon_coord_w3nt_dict = {}
    protein_seq_dict = {}
    cds_seq_dict = {}
    cds_exon_coord_list_dict = {}
    for tx_id in tx_list:
        chrom = chrom_dict[tx_id]
        strand = strand_dict[tx_id]
        exon_coord_list = tx_exon_coord_dict[tx_id]
        tx_seq_w3nt, exon_coord_list_w3nt = get_seq_from_exon_list(chrom, strand, exon_coord_list, add_tail=3)
        if tx_seq_w3nt == '':
            if verbose == 1:
                print('%s does not have valid exon coordinates, excluded.'%(tx_id))
            continue
        protein_seq, cds_seq, cds_exon_coord_list, cds_coord_in_tx = translate_tx(tx_seq_w3nt, strand, exon_coord_list_w3nt)
        protein_seq = protein_seq.rstrip('*')
        if protein_seq == '':
            if verbose == 1:
                print('%s does not have in-frame ORF, excluded.'%(tx_id))
            continue
        if protein_seq.find('*') >= 0: #ambiguous PTC exists in the sequences
            if verbose == 1:
                print('%s: protein sequence has PTC, excluded'%(tx_id))
            continue
        if protein_seq.find('X') >= 0: #ambiguous PTC exists in the sequences
            if verbose == 1:
                print('%s: protein sequence has unknown amino acid X, excluded'%(tx_id))
            continue
        tx_seq_w3nt_dict[tx_id] = tx_seq_w3nt
        tx_exon_coord_w3nt_dict[tx_id] = exon_coord_list_w3nt
        protein_seq_dict[tx_id] = protein_seq
        cds_seq_dict[tx_id] = cds_seq
        cds_exon_coord_list_dict[tx_id] = cds_exon_coord_list
        tx_list_wCDS.append(tx_id)

    outf_gtf = open(output_gtf, 'w')
    if output_transcript:
        outf_transcript = open(output_transcript, 'w')
    else:
        outf_transcript = None

    if output_cds:
        outf_cds = open(output_cds, 'w')
    else:
        outf_cds = None

    if output_protein:
        outf_protein = open(output_protein, 'w')
    else:
        outf_protein = None

    for tx_id in tx_list_wCDS:
        chrom = chrom_dict[tx_id]
        strand = strand_dict[tx_id]
        exon_coord_list = tx_exon_coord_dict.get(tx_id, [])
        exon_coord_list_w3nt = tx_exon_coord_w3nt_dict.get(tx_id, [])
        cds_exon_coord_list = cds_exon_coord_list_dict.get(tx_id, [])
        tx_seq_w3nt = tx_seq_w3nt_dict[tx_id]
        cds_seq = cds_seq_dict[tx_id]
        _NMD_result = check_NMD_by_seq(tx_seq_w3nt, cds_seq)
        if _NMD_result:
            if keep_nmd == 'True':
                tx_id_wNMD = '%s'%tx_id
            elif keep_nmd == 'Add':
                tx_id_wNMD = '%s NMD'%tx_id
            elif keep_nmd == 'False':
                if verbose == 1:
                    print('%s is predicted to have NMD, excluded.'%(tx_id))
                continue
        else:
            if keep_nmd == 'True':
                if verbose == 1:
                    print('Protein coding transcript %s is excluded when keep_nmd is set to True.'%(tx_id))
                continue
            elif keep_nmd == 'Add':
                tx_id_wNMD = '%s'%tx_id
            elif keep_nmd == 'False':
                tx_id_wNMD = '%s'%tx_id
        protein_seq = protein_seq_dict[tx_id]
        if outf_protein:
            outf_protein.write('>%s\n%s\n'%(tx_id_wNMD, protein_seq))
        if outf_transcript:
            outf_transcript.write('>%s\n%s\n'%(tx_id_wNMD, tx_seq_w3nt))
        if outf_cds:
            outf_cds.write('>%s\n%s\n'%(tx_id_wNMD, cds_seq))
        if outf_gtf:
            gtf_line = tx_gtf_dict.get(tx_id, '')
            if gtf_line != '':
                outf_gtf.write(gtf_line)
            d = filter_attributes(tx_attribute_dict.get(tx_id, {}), attribute_list)
            attribute = deparse_attributes(d)
            exon_coord_list_adjusted_3nt = compare_exon_coords_vs_cds_coords(exon_coord_list, cds_exon_coord_list, strand)
            for exon_coord in exon_coord_list_adjusted_3nt:
                gtf_line = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n'%(chrom, 'ProteoSeq', 'exon', exon_coord[0], exon_coord[1], '.', strand, '.', attribute)
                outf_gtf.write(gtf_line)
            for cds_coord in cds_exon_coord_list:
                gtf_line = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n'%(chrom, 'ProteoSeq', 'CDS', cds_coord[0], cds_coord[1], '.', strand, '.', attribute)
                outf_gtf.write(gtf_line)

    if outf_gtf:
        outf_gtf.close()
    if outf_transcript:
        outf_transcript.close()
    if outf_cds:
        outf_cds.close()
    if outf_protein:
        outf_protein.close()

