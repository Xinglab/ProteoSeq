import os,sys,argparse
import random

def reverse_seq(seq):
    return seq[::-1]

def shuffle_seq(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return ''.join(seq_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reverse or/and shuffle protein sequences.')
    parser.add_argument('-i', '--input_file', help='Input FASTA. Required', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output FASTA. Required', type=str, required=True)
    parser.add_argument('--prefix', help='Rename output protein sequences as >prefix_N. If none, keep original name', type=str)
    parser.add_argument('--reverse', help='Set to reverse each sequence.', action='store_true')
    parser.add_argument('--shuffle', help='Set to shuffle each sequence.', action='store_true')
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    prefix = ''
    if args.prefix:
        prefix = args.prefix
    reverse = False
    if args.reverse:
        reverse = True
    shuffle = False
    if args.shuffle:
        shuffle = True

    protein_seq_dict = {}
    n = 0
    for line in open(input_file, 'r'):
        if line[0] == '>':
            seq_id = line.rstrip('\n')[1:]
            n += 1
        else:
            protein_seq_dict[(n, seq_id)] = protein_seq_dict.get((n, seq_id), '') + line.rstrip('\n')

    output = open(output_file, 'w')
    for (n, seq_id) in sorted(protein_seq_dict):
        seq = protein_seq_dict[(n, seq_id)]
        if reverse:
            seq = reverse_seq(seq)
        if shuffle:
            seq = shuffle_seq(seq)
        if prefix:
            output.write('>%s_%d\n%s\n'%(prefix, n, seq))
        else:
            output.write('>%s\n%s\n'%(seq_id, seq))

    output.close()



