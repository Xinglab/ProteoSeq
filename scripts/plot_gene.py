import os,re,sys,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Rectangle

exon_height = .7
font_size = 9
arrow_length = .3

def parse_attributes(string):
    d = {}
    for i in string.split(';'):
        if i == '':
            continue
        arr = i.strip().split(' ')
        key = arr[0]
        value = ' '.join(arr[1:]).strip('"')
        d.setdefault(key, []).append(value)
    return d

def get_exon_genomic_coords(strand, tx_start, tx_end, exon_start_s0_list, exon_len_list):
    exon_coord_list_genomic = []
    for i in range(len(exon_start_s0_list)):
        exon_start_s0 = exon_start_s0_list[i]
        exon_len = exon_len_list[i]
        exon_start_genomic = tx_start + exon_start_s0
        exon_end_genomic = exon_start_genomic + exon_len - 1
        exon_coord_list_genomic.append((exon_start_genomic, exon_end_genomic))
    return exon_coord_list_genomic

def get_exon_structures(bed_list):
    exon_dict = {}
    chrom_dict = {}
    strand_dict = {}
    for arr in bed_list:
        chrom, tx_id, strand,  = arr[0], arr[3], arr[5]
        tx_start_s0, tx_end = int(arr[1]), int(arr[2])
        tx_start = tx_start_s0 + 1
        exon_len_list = [int(_) for _ in arr[10].rstrip(',').split(',')]
        exon_start_s0_list = [int(_) for _ in arr[11].rstrip(',').split(',')]
        strand_dict[tx_id] = strand
        chrom_dict[tx_id] = chrom
        exon_coord_list_genomic = get_exon_genomic_coords(strand, tx_start, tx_end, exon_start_s0_list, exon_len_list)
        exon_dict[tx_id] = exon_coord_list_genomic
    return exon_dict, chrom_dict, strand_dict

def adjust_pos(exon_dict, min_pos):
    final_pos_dict = {}
    corrected_exon_dict = {}
    corrected_intron_dict = {}
    for tx_id, exon_coord_list in exon_dict.items():
        if len(exon_coord_list) == 0:
            corrected_exon_dict[tx_id] = []
            corrected_intron_dict[tx_id] = []
            final_pos_dict[tx_id] = 0
            continue
        exon_n = len(exon_coord_list)
        ### adjust initial current_pos
        initial_pos = exon_coord_list[0][0] - min_pos
        current_pos = initial_pos
        for i in range(exon_n):
            exon_coord = exon_coord_list[i]
            exon_len = exon_coord[1] - exon_coord[0] + 1
            current_exon_pos = (current_pos, current_pos+exon_len)
            corrected_exon_dict.setdefault(tx_id, []).append(current_exon_pos)
            current_pos += exon_len
            if (exon_n >= 2) and ( i < (exon_n-1)):
                exon_coord_next = exon_coord_list[i+1]
                intron_len = exon_coord_next[0] - exon_coord[1] - 1
                current_intron_pos = (current_pos, current_pos+intron_len)
                corrected_intron_dict.setdefault(tx_id, []).append(current_intron_pos)
                current_pos = current_pos + intron_len
        final_pos_dict[tx_id] = current_pos
    return final_pos_dict,corrected_exon_dict,corrected_intron_dict

def get_min_max_pos(exon_dict):
    min_pos, max_pos = 10000000000, 0
    for tx_id, exon_coord_list in exon_dict.items():
        if len(exon_coord_list) == 0:
            continue
        min_pos = min(min_pos, min([_[0] for _ in exon_coord_list]))
        max_pos = max(max_pos, max([_[1] for _ in exon_coord_list]))
    return min_pos, max_pos

def select_exons(exon_coord_list, select_start, select_end):
    selected_exon_coord_list = []
    for exon_coord in exon_coord_list:
        if exon_coord[1] < select_start:
            continue
        if exon_coord[0] > select_end:
            continue
        selected_exon_coord_list.append(exon_coord)
    return selected_exon_coord_list

def select_color(tx_id):
    if tx_id.startswith('ENST'):
        return '#003E7F'
    return '#FF5F42'

def select_hatch(tx_id):
    return None

def draw_intron(ax, y_pos, exon_coord_list, strand, color='r'):
    if len(exon_coord_list) < 2:
        return False
    intron_left = exon_coord_list[0][1]
    intron_right = exon_coord_list[-1][0]
    if strand == '+':
        ax.hlines(y_pos, intron_left, intron_right, color=color,lw=2)
    else:
        ax.hlines(y_pos, -intron_right, -intron_left, color=color,lw=2)

def draw_exon(ax, y_pos, exon_coord_list, strand, color='r', hatch=None):
    if len(exon_coord_list) == 0:
        return False
    for (exon_left, exon_right) in exon_coord_list:
        exon_length = exon_right - exon_left
        if strand == '+':
            ax.add_patch(Rectangle((exon_left, y_pos-.5*exon_height), exon_length, exon_height, fill=True, fc=color, lw=0, hatch=hatch))
        else:
            ax.add_patch(Rectangle((-exon_right, y_pos-.5*exon_height), exon_length, exon_height, fill=True, fc=color, lw=0, hatch=hatch))

def plot_tx_structure(ax, tx_id, exon_coord_list, y_pos, strand, color, hatch=None):
    draw_intron(ax, y_pos, exon_coord_list, strand, color)
    draw_exon(ax, y_pos, exon_coord_list, strand, color, hatch=hatch)
    
def plot_structures(tx_list, exon_coord_dict, final_pos_dict, color_list=[], tx_id_label_list=[], output_file=None):
    num_trans = len(tx_list)
    fig, ax = plt.subplots(1, 1, figsize=(10, .5*num_trans))
    ymin = .5
    ymax = num_trans+1
    ax.set_ylim(ymin, ymax)
    
    max_final_pos = 0
    for tx_idx in range(num_trans):
        y_pos = num_trans - tx_idx
        tx_id = tx_list[tx_idx]
        strand = tx_strand_dict[tx_id]
        final_pos = final_pos_dict[tx_id]
        max_final_pos = max(max_final_pos, final_pos)
        if len(color_list) > 0 and tx_idx < len(color_list):
            color = color_list[tx_idx]
        else:
            color = select_color(tx_id)
        hatch = select_hatch(tx_id)
        exon_coord_list = exon_coord_dict[tx_id]
        plot_tx_structure(ax, tx_id, exon_coord_list, y_pos, strand, color, hatch=hatch)
    xmin, xmax = max_final_pos*.01, max_final_pos*1.01
    if strand == '+':
        ax.set_xlim(-xmin, xmax)
    else:
        ax.set_xlim(-xmax, xmin)
    for tx_idx in range(num_trans):
        y_pos = num_trans - tx_idx
        tx_id = tx_list[tx_idx]
        strand = tx_strand_dict[tx_id]
        if tx_id_label_list and len(tx_id_label_list) >= num_trans:
            tx_id_label = tx_id_label_list[tx_idx]
        else:
            tx_id_label = tx_id
        if strand == '+':
            ax.text(0, y_pos+exon_height/2, tx_id_label, va='bottom', ha='left', fontsize=font_size+1)
            ax.text(-xmin, y_pos, '5\'', va='center', ha='right', fontsize=font_size+2)
        else:
            ax.text(-max_final_pos, y_pos+exon_height/2, tx_id_label, va='bottom', ha='left', fontsize=font_size+1)
            ax.text(-xmax, y_pos, '5\'', va='center', ha='right', fontsize=font_size+2)
        
    
    ax.text(.5, 1, gene_name, ha='center', va='top', fontsize=font_size+2, transform=ax.transAxes)
    ax.axis('off')
    
    fig.subplots_adjust(top=1, bottom=0, left=.03, right=.985)
    if output_file:
        fig.savefig(output_file) #, dpi=300
    plt.close()
    plt.clf()

def get_tx_id_label(tx_id):
    tx_id_label_list = []
    if tx_canonical_dict.get(tx_id, 0) > 0:
        tx_id_label_list.append('Canonical')
    tx_type = tx_type_dict.get(tx_id, None)
    if tx_type == 'protein_coding':
        tx_id_label_list.append('Protein_coding')
    if tx_basic_dict.get(tx_id, None):
        tx_id_label_list.append('Basic')
    if len(tx_id_label_list) > 0:
        return '%s (%s)'%(tx_id, '; '.join(tx_id_label_list))
    return tx_id
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot transcripts structures')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file of transcripts and exon coordinates. Required', type=str, required=True)
    parser.add_argument('-g', '--gene_name_list', help='Plot for every given gene', nargs='+', type=str, required=True)
    parser.add_argument('-o', '--output_folder', help='Output folder to plot figures', type=str, required=True)
    args = parser.parse_args()

    input_gtf = args.input_gtf
    gene_name_list = args.gene_name_list
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    tx_canonical_dict = {}
    tx_basic_dict = {}
    tx_type_dict = {}
    tx_basic_dict = {}
    gene_id_to_name_dict = {}
    detected_gene_id_to_tx_dict = {}
    detected_tx_to_gene_id_dict = {}
    tx_coord_dict = {}
    tx_exon_coord_dict = {}
    tx_chrom_dict = {}
    tx_strand_dict = {}
    for line in open(input_gtf, 'r'):
        if line[0] == '#':
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            gene_id = d.get('gene_id', [''])[0]
            if gene_id == '':
                continue
            gene_name = d.get('gene_name', [''])[0]
            if gene_name not in gene_name_list:
                continue
            detected_tx_to_gene_id_dict[tx_id] = gene_id
            detected_gene_id_to_tx_dict.setdefault(gene_id, []).append(tx_id)
            gene_id_to_name_dict[gene_id] = gene_name
            tx_coord = (int(arr[3]), int(arr[4]))
            tx_coord_dict[tx_id] = tx_coord
            tx_chrom_dict[tx_id] = arr[0]
            tx_strand_dict[tx_id] = arr[6]
            tags = d.get('tag', [])
            if 'Ensembl_canonical' in tags:
                tx_canonical_dict[tx_id] = 1
            if 'basic' in tags:
                tx_basic_dict[tx_id] = 'basic'
            tx_type = d.get('transcript_type', [''])[0]
            if tx_type:
                tx_type_dict[tx_id] = tx_type
        elif arr[2] == 'exon':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            exon_coord = (int(arr[3]), int(arr[4]))
            tx_exon_coord_dict.setdefault(tx_id, []).append(exon_coord)

    for gene_id in detected_gene_id_to_tx_dict.keys():
        gene_name = gene_id_to_name_dict.get(gene_id, gene_id)
        detected_tx_list = detected_gene_id_to_tx_dict[gene_id]
        tx_id_label_list = []
        for tx_id in detected_tx_list:
            tx_id_label = get_tx_id_label(tx_id)
            tx_id_label_list.append(tx_id_label)
        num_trans = len(detected_tx_list)
        exon_dict = {tx_id:tx_exon_coord_dict[tx_id] for tx_id in detected_tx_list}
        min_pos, max_pos = get_min_max_pos(exon_dict)
        final_pos_dict,corrected_exon_dict,corrected_intron_dict = adjust_pos(exon_dict, min_pos)
        plot_structures(detected_tx_list, corrected_exon_dict, final_pos_dict,
                        tx_id_label_list=tx_id_label_list,
                        output_file='%s/%s_structure.png'%(output_folder, gene_name))

