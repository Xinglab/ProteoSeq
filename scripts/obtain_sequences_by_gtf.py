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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Obtain transcript sequences with annotated coordinates in a GTF.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file with coordinates. Required', type=str, required=True)
    parser.add_argument('-g', '--genome_fasta', help='Input genome fasta file. Required', type=str, required=True)
    parser.add_argument('-o', '--output_transcript', help='Output fasta file of transcript sequences. Required', type=str, required=True)
    parser.add_argument('-v', '--verbose', help='Set to zero to silent stdout skipped lines, default=0', type=int, choices=[0, 1], default=0)
    args = parser.parse_args()
    
    input_gtf = args.input_gtf
    output_transcript = args.output_transcript
    genome_fasta = args.genome_fasta
    verbose = args.verbose

    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    tx_exon_coord_dict = {}
    chrom_dict = {}
    strand_dict = {}
    for line in open(input_gtf, 'r'):
        if line.startswith('#'):
            continue
        arr = line.rstrip('\n').split('\t')
        if len(arr) < 2:
            #print(line)
            continue
        if arr[2] == 'exon':
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
    tx_list = sorted(tx_exon_coord_dict)
    
    tx_seq_dict = {}
    for tx_id in tx_list:
        chrom = chrom_dict[tx_id]
        strand = strand_dict[tx_id]
        exon_coord_list = tx_exon_coord_dict.get(tx_id, [])
        
        if len(exon_coord_list) == 0:
            if verbose == 1:
                print('%s: protein sequence has unknown amino acid X, excluded'%(tx_id))
            continue
        tx_seq = get_seq_from_exon_list(chrom, strand, exon_coord_list)
        tx_seq_dict[tx_id] = tx_seq

    outf_transcript = open(output_transcript, 'w')
    for tx_id,tx_seq in tx_seq_dict.items():
         _ = outf_transcript.write('>%s\n%s\n'%(tx_id, tx_seq))

    outf_transcript.close()


