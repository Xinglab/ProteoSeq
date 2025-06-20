"""
Script Name: select_species_in_uniprot.py
Description: Select sequences in certain species with flag OX=TAXA_NUM in the sequence IDs.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,re,sys,argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select sequences in certain species with flag OX=TAXA_NUM in the sequence IDs.')
    parser.add_argument('-i', '--input_fasta', help='Input FASTA file of sequences. Required', type=str, required=True)
    parser.add_argument('-o', '--output_fasta', help='Output FASTA', type=str, required=True)
    parser.add_argument('--species', help='Taxanomic numbers of species retained in the FASTA file, comma splitted of multiple spieces. default=9606', type=str, default='9606')

    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_fasta = args.output_fasta
    species_list = set(args.species.split(','))
    
    seq_dict = {}
    for line in open(input_fasta, 'r'):
        if line.startswith('>'):
            seq_id = line.rstrip('\n')
        else:
            seq_dict[seq_id] = seq_dict.get(seq_id, '') + line.rstrip('\n')


    output = open(output_fasta, 'w')
    for seq_id, seq in seq_dict.items():
        species = re.findall('OX=(\d+)', seq_id)
        if len(species) == 0:
            continue
        if species[0] in species_list:
            output.write('%s\n%s\n'%(seq_id, seq))

    output.close()
