"""
Script Name: translate_by_gtf.py
Description: Translate protein sequences with annotated CDS coordinates in a GTF.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,sys,argparse
import re
from Bio.Seq import Seq
from Bio import SeqIO

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
        seq = temp_seq[:-1] + temp_seq[-1].lower()
        final_seq += seq
    temp_final_seq = final_seq[:-1] + final_seq[-1].upper()
    final_seq = temp_final_seq
    return final_seq

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

def lower_ss(seq, pos_s0_list):
    if len(pos_s0_list)==0:
        return seq
    for pos_s0 in pos_s0_list:
        if 0 <= pos_s0 < len(seq):
            temp_seq = seq[:pos_s0] + seq[pos_s0].lower() + seq[pos_s0+1:]
            seq = temp_seq
    return seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Translate protein sequences with annotated CDS coordinates in a GTF.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file with CDS coordinates. Required', type=str, required=True)
    parser.add_argument('-g', '--genome_fasta', help='Input genome fasta file. Required', type=str, required=True)
    parser.add_argument('-o', '--output_protein', help='Output fasta file of protein sequences. Required', type=str, required=True)
    parser.add_argument('--output_transcript', help='Output additional fasta file of transcript sequences', type=str)
    parser.add_argument('--output_cds', help='Output additional fasta file of CDS sequences', type=str)
    parser.add_argument('--keep_ptc', help='Keep sequences with PTC (truncated at the first stop codon), default=False ', type=str, default='False', choices=['True', 'False'])
    parser.add_argument('--not_triplet', help='Force translation of ORFs with length not multiple of three, default=False ', type=str, default='False', choices=['True', 'False'])
    parser.add_argument('-v', '--verbose', help='Set to zero to silent stdout skipped lines, default=0', type=int, choices=[0, 1], default=0)
    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_protein = args.output_protein
    genome_fasta = args.genome_fasta

    if args.output_transcript:
        output_transcript = args.output_transcript
    else:
        output_transcript = None

    if args.output_cds:
        output_cds = args.output_cds
    else:
        output_cds = None

    if args.keep_ptc == 'True':
        keep_ptc = True
    else:
        keep_ptc = False
        
    if args.not_triplet == 'True':
        not_triplet = True
    else:
        not_triplet = False
    
    verbose = args.verbose
    
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    
    tx_exon_coord_dict = {}
    tx_cds_coord_dict = {}
    chrom_dict = {}
    strand_dict = {}
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if len(arr) < 2:
            #print(line)
            continue
        if arr[2] == 'CDS':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('Excluded line in input GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            tx_cds_coord_dict.setdefault(tx_id, []).append((int(arr[3]), int(arr[4])))
            chrom_dict[tx_id] = arr[0]
            strand_dict[tx_id] = arr[6]
        elif arr[2] == 'exon':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                if verbose == 1:
                    print('Excluded line in input GTF (no transcript id): %s'%(line.rstrip('\n')))
                continue
            tx_exon_coord_dict.setdefault(tx_id, []).append((int(arr[3]), int(arr[4])))
            chrom_dict[tx_id] = arr[0]
            strand_dict[tx_id] = arr[6]

    tx_exon_coord_dict = {k:sorted(v) for k,v in tx_exon_coord_dict.items()}
    tx_cds_coord_dict = {k:sorted(v) for k,v in tx_cds_coord_dict.items()}



    tx_list_wCDS = sorted(tx_cds_coord_dict)
    
    tx_seq_dict = {}
    cds_seq_dict = {}
    protein_seq_dict = {}
    for tx_id in tx_list_wCDS:
        chrom = chrom_dict[tx_id]
        strand = strand_dict[tx_id]
        cds_exon_coord_list = tx_cds_coord_dict.get(tx_id, [])
        exon_coord_list = tx_exon_coord_dict.get(tx_id, [])
        if len(cds_exon_coord_list) == 0: #No CDS is found in annotation
            continue
        cds_seq = get_seq_from_exon_list(chrom, strand, cds_exon_coord_list)
        protein_seq = Seq(cds_seq).translate().rstrip('*')
        if not protein_seq:
            if verbose == 1:
                print('%s: invalid CDS coordinates, excluded'%(tx_id))
            continue
        if len(cds_seq) % 3 > 0: #not multiple of three
            if not not_triplet:
                if verbose == 1:
                    print('%s: CDS length is not multiple of three, excluded'%(tx_id))
                continue
        if protein_seq.find('*') >= 0: #ambiguous PTC exists in the sequences
            if keep_ptc:
                protein_seq = protein_seq[:protein_seq.find('*')] #truncate at the first stop codon
                if not protein_seq:
                    continue
            else:
                if verbose == 1:
                    print('%s: protein sequence has PTC, excluded'%(tx_id))
                continue
        if protein_seq.find('X') >= 0: #ambiguous PTC exists in the sequences
            if keep_ptc:
                protein_seq = protein_seq[:protein_seq.find('*')] #truncate at the first stop codon
                if not protein_seq:
                    continue
            else:
                if verbose == 1:
                    print('%s: protein sequence has unknown amino acid X, excluded'%(tx_id))
                continue
        ss_pos_s0_in_cds = [i.start() for i in re.finditer("[a-z]", cds_seq)]
        ss_pos_s0_in_aa = [_//3 for _ in ss_pos_s0_in_cds]
        protein_seq_with_ss = lower_ss(protein_seq, ss_pos_s0_in_aa)
        if protein_seq_with_ss[-1]=='*': #remove final stop codon-derived AA
            protein_seq_with_ss = protein_seq_with_ss[0:-1]
        protein_seq_dict[tx_id] = protein_seq_with_ss
        cds_seq_dict[tx_id] = cds_seq
        if len(exon_coord_list) > 0:
            tx_seq = get_seq_from_exon_list(chrom, strand, exon_coord_list)
            tx_seq_dict[tx_id] = tx_seq


    outf_protein = open(output_protein, 'w')
    
    if output_transcript:
        outf_transcript = open(output_transcript, 'w')
    else:
        outf_transcript = None

    if output_cds:
        outf_CDS = open(output_cds, 'w')
    else:
        outf_CDS = None

    for tx_id,protein_seq in protein_seq_dict.items():
        _ = outf_protein.write('>%s\n%s\n'%(tx_id, protein_seq))
        if outf_CDS:
            cds_seq = cds_seq_dict[tx_id]
            _ = outf_CDS.write('>%s\n%s\n'%(tx_id, cds_seq))
        if outf_transcript:
            tx_seq = tx_seq_dict[tx_id]
            _ = outf_transcript.write('>%s\n%s\n'%(tx_id, tx_seq))

    outf_protein.close()
    
    if outf_transcript:
        outf_transcript.close()

    if outf_CDS:
        outf_CDS.close()


