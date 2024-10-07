import os,re,sys,argparse
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter out protein sequences with premature stop codon (\"*\") and instandard amino acids (e.g., \"X\" in the sequences).')
    parser.add_argument('-i', '--input_fasta', help='Input fasta file of protein squences. Required', type=str, required=True)
    parser.add_argument('-o', '--output_fasta', help='Output fasta file of protein sequences w/o PTC. Required', type=str, required=True)
    parser.add_argument('--output_ptc', help='Output additional fasta file of protein sequences containing PTC', type=str)
    parser.add_argument('--output_nmd', help='Output additional fasta file of protein sequences containing NMD. Setting this on, the script will look for sequence having \"NMD\" in their IDs.', type=str)
    parser.add_argument('--output_ambiguous', help='Output additional fasta file of protein sequences containing instandard amino acids such as \"X\".', type=str)
    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_fasta = args.output_fasta

    if args.output_ptc:
        output_ptc = args.output_ptc
    else:
        output_ptc = None

    if args.output_nmd:
        output_nmd = args.output_nmd
    else:
        output_nmd = None

    if args.output_ambiguous:
        output_ambiguous = args.output_ambiguous
    else:
        output_ambiguous = None

    pc_seq_dict = {}
    nmd_seq_dict = {}
    ptc_seq_dict = {}
    ambiguous_seq_dict = {}
    nmd_flag = False
    for line in open(input_fasta, 'r'):
        if line[0] == '>':
            nmd_flag = False
            if re.findall('NMD', line) and (output_nmd != None):
                nmd_flag = True
            seq_id = line.rstrip('\n')[1:]
        else:
            seq = line.rstrip('\n').upper()
            if nmd_flag:
                nmd_seq_dict[seq_id] = nmd_seq_dict.get(seq_id, '') + seq
            elif (re.findall('\*', seq)):
                ptc_seq_dict[seq_id] = ptc_seq_dict.get(seq_id, '') + seq
            elif re.findall('X', seq):
                ambiguous_seq_dict[seq_id] = ambiguous_seq_dict.get(seq_id, '') + seq
            else:
                pc_seq_dict[seq_id] = pc_seq_dict.get(seq_id, '') + seq

    output = open(output_fasta,'w')
    for seq_id, seq in pc_seq_dict.items():
        output.write('>%s\n%s\n' % (seq_id, seq))
    output.close()

    if output_nmd:
        output = open(output_nmd,'w')
        for seq_id, seq in nmd_seq_dict.items():
            output.write('>%s\n%s\n' % (seq_id, seq))
        output.close()
    else:
        for seq_id, seq in nmd_seq_dict.items():
            print('%s is excluded due to NMD in the sequence'%(seq_id))

    if output_ptc:
        output = open(output_ptc,'w')
        for seq_id, seq in ptc_seq_dict.items():
            output.write('>%s\n%s\n' % (seq_id, seq))
        output.close()
    else:
        for seq_id, seq in ptc_seq_dict.items():
            print('%s is excluded due to PTC in the sequence'%(seq_id))

    if output_ambiguous:
        output = open(output_ambiguous,'w')
        for seq_id, seq in ambiguous_seq_dict.items():
            output.write('>%s\n%s\n' % (seq_id, seq))
        output.close()
    else:
        for seq_id, seq in ambiguous_seq_dict.items():
            print('%s is excluded due to ambiguous amino acids in the sequence'%(seq_id))


