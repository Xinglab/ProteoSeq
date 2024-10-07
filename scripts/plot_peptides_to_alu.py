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


def select_color(tx_id):
    if tx_id.startswith('ENST'):
        return 'blue'
    return 'red'

def draw_intron(ax, y_pos, exon_coord_list, strand, color='r'):
    if len(exon_coord_list) < 2:
        return False
    intron_left = exon_coord_list[0][1]
    intron_right = exon_coord_list[-1][0]
    if strand == '+':
        ax.hlines(y_pos, intron_left, intron_right, color=color,lw=2)
    else:
        ax.hlines(y_pos, -intron_right, -intron_left, color=color,lw=2)

def draw_exon(ax, y_pos, exon_coord_list, final_pos, strand, color='r', hatch=None):
    if len(exon_coord_list) == 0:
        return False
    for (exon_left, exon_right) in exon_coord_list:
        exon_length = exon_right - exon_left
        if exon_length/final_pos < 0.001:
            exon_length = final_pos*0.001
        if strand == '+':
            ax.add_patch(Rectangle((exon_left, y_pos-.5*exon_height), exon_length, exon_height, fill=True, fc=color, lw=0, hatch=hatch))
        else:
            ax.add_patch(Rectangle((-exon_right, y_pos-.5*exon_height), exon_length, exon_height, fill=True, fc=color, lw=0, hatch=hatch))

def plot_peptide_arrows(ax, y_pos, peptide_start, peptide_end, strand, max_final_pos=0, text=None):
    if strand == '+':
        ax.annotate("", xy=(peptide_start, y_pos), xytext=(peptide_start, y_pos-arrow_length), arrowprops=dict(arrowstyle="->"))
        ax.annotate("", xy=(peptide_end, y_pos), xytext=(peptide_end, y_pos-arrow_length), arrowprops=dict(arrowstyle="->"))
        if text:
            ax.text(peptide_start-max_final_pos*.008, y_pos-arrow_length, text, ha='right', va='bottom', fontsize=font_size-2)
    else:
        ax.annotate("", xy=(-peptide_start, y_pos), xytext=(-peptide_start, y_pos-arrow_length), arrowprops=dict(arrowstyle="->"))
        ax.annotate("", xy=(-peptide_end, y_pos), xytext=(-peptide_end, y_pos-arrow_length), arrowprops=dict(arrowstyle="->"))
        if text:
            ax.text(-peptide_end-max_final_pos*.008, y_pos-arrow_length, text, ha='right', va='bottom', fontsize=font_size-2)


def plot_tx_structure(ax, tx_id, exon_coord_list, final_pos, y_pos, strand, color, hatch=None, text=None):
    draw_intron(ax, y_pos, exon_coord_list, strand, color)
    draw_exon(ax, y_pos, exon_coord_list, final_pos, strand, color, hatch=hatch)
    if text:
        (exon_left, exon_right) = exon_coord_list[0]
        if strand == '+':
            ax.text(exon_left, y_pos+.5*exon_height, text, va='bottom', ha='left', fontsize=font_size-2, rotation=90)
        else:
            ax.text(-exon_right, y_pos+.5*exon_height, text, va='bottom', ha='left', fontsize=font_size-2, rotation=90)


def plot_structures(tx_list, te_list, exon_coord_dict, final_pos_dict, peptide_count_dict=None, output_file=None):
    num_trans = len(tx_list)
    num_te = len(te_list)
    num_peptide = 0
    peptide_coord_list = []
    if peptide_count_dict:
        peptide_coord_list = sorted(peptide_count_dict, key=lambda k:peptide_count_dict[k], reverse=True)
        num_peptide = len(peptide_coord_list)
    fig, ax = plt.subplots(1, 1, figsize=(10, .5*(num_trans + arrow_length*num_peptide + 2)))
    ymin = .5-arrow_length*num_peptide
    ymax = num_trans + 3
    ax.set_ylim(ymin, ymax)
    max_final_pos = 0
    for tx_idx in range(num_trans):
        y_pos = tx_idx + 1
        tx_id = tx_list[tx_idx]
        strand = tx_strand_dict[tx_id]
        final_pos = final_pos_dict[tx_id]
        max_final_pos = max(max_final_pos, final_pos)
        color = select_color(tx_id)
        exon_coord_list = exon_coord_dict[tx_id]
        plot_tx_structure(ax, tx_id, exon_coord_list, final_pos, y_pos, strand, color)
        ax.text(0, (y_pos+exon_height*.5-ymin)/(ymax-ymin), tx_id, va='bottom', ha='left', fontsize=font_size+1, transform=ax.transAxes)
        ax.text(0, (y_pos-ymin)/(ymax-ymin), '5\'', va='center', ha='right', fontsize=font_size+2, transform=ax.transAxes)
    if strand == '+':
        ax.set_xlim(0, max_final_pos)
    else:
        ax.set_xlim(-max_final_pos, 0)
    for te_idx in range(num_te):
        y_pos = num_trans + 1
        te_id = te_list[te_idx]
        final_pos = final_pos_dict[te_id]
        max_final_pos = max(max_final_pos, final_pos)
        exon_coord_list = corrected_exon_dict[te_id]
        plot_tx_structure(ax, te_id, exon_coord_list, final_pos, y_pos, strand, 'gray', text=te_id.split(' ')[0])
    if num_peptide > 0:
        for peptide_idx in range(num_peptide):
            peptide_coord = peptide_coord_list[peptide_idx]
            peptide_count = peptide_count_dict[peptide_coord]
            plot_peptide_arrows(ax, .5-arrow_length*peptide_idx, peptide_coord[0], peptide_coord[1], strand, max_final_pos, text=str(peptide_count))
    ax.text(.5, 1, gene_name, ha='center', va='top', fontsize=font_size+2, transform=ax.transAxes)
    ax.axis('off')
    fig.subplots_adjust(top=1, bottom=0, left=.03, right=.985)
    if output_file:
        fig.savefig(output_file, dpi=300)
    plt.close()
    plt.clf()


def check_overlap(coord1, coord2):
    sorted_coord = sorted([coord1[0], coord1[1], coord2[0], coord2[1]])
    if sorted_coord == [coord1[0], coord1[1], coord2[0], coord2[1]]:
        return False
    if sorted_coord == [coord2[0], coord2[1], coord1[0], coord1[1]]:
        return False
    return True


def check_te_with_tx_exons(te_coord, tx_exon_coord_list):
    for exon_coord in tx_exon_coord_list:
        if check_overlap(te_coord, exon_coord):
            return True
    return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot peptides to the transcripts structures with Alu exons')
    parser.add_argument('-p', '--peptide_table', help='Input peptide table. Required', type=str, required=True)
    parser.add_argument('-r', '--ref_gtf', help='Input GTF file of protein references. Required', type=str, required=True)
    parser.add_argument('-o', '--output_folder', help='Output folder to plot figures', type=str, required=True)
    parser.add_argument('--tx_bed', help='Bed of transcripts with peptides', type=str, required=True)
    parser.add_argument('--alu_bed', help='Bed of transcripts with peptides', type=str, required=True)
    args = parser.parse_args()

    peptide_table = args.peptide_table
    ref_gtf = args.ref_gtf
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    tx_bed = args.tx_bed
    alu_bed = args.alu_bed

    peptide_count_dict = {}
    header_dict = {}
    for line in open(peptide_table, 'r'):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        peptide_seq = arr[0]
        psm_count = int(arr[header_dict['total_count']])
        peptide_count_dict[peptide_seq] = peptide_count_dict.get(peptide_seq, 0) + psm_count

    tx_matched_te_bed_dict = {}
    for line in open(tx_bed, 'r'):
        arr = line.rstrip('\n').split('\t')
        tx_id = arr[3].split('|')[0]
        te_bed = tuple(arr[6:12])
        tx_matched_te_bed_dict.setdefault(tx_id, set()).add(te_bed)

    tx_peptides_bed_dict = {}
    for line in open(alu_bed, 'r'):
        arr = line.rstrip('\n').split('\t')
        arr2 = arr[3].split('_')
        peptide_seq = arr2[0]
        tx_id_list = '_'.join(arr2[1:]).split(';')
        peptide_coord = (arr[0], arr[5], int(arr[1]), int(arr[2]))
        #te_bed = tuple(arr[6:12])
        for tx_id in tx_id_list:
            tx_peptides_bed_dict.setdefault(tx_id, {})[peptide_coord] = peptide_seq
        
    gene_id_to_name_dict = {}
    detected_tx_to_gene_id_dict = {}
    tx_exon_coord_dict = {}
    tx_strand_dict = {}
    for line in open(ref_gtf, 'r'):
        if line[0] == '#':
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            gene_id = d.get('gene_id', [''])[0]
            gene_name = d.get('gene_name', [''])[0]
            if gene_id != '':
                detected_tx_to_gene_id_dict[tx_id] = gene_id
                if gene_name != '':
                    gene_id_to_name_dict[gene_id] = gene_name
            tx_strand_dict[tx_id] = arr[6]
        elif arr[2] == 'exon':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            exon_coord = (int(arr[3]), int(arr[4]))
            tx_exon_coord_dict.setdefault(tx_id, []).append(exon_coord)

    gene_to_tx_list_with_alu_peptides_dict = {}
    for tx_id in sorted(tx_peptides_bed_dict.keys()):
        bed_list_te = tx_matched_te_bed_dict.get(tx_id, None)
        if not bed_list_te:
            continue
        gene_id = detected_tx_to_gene_id_dict.get(tx_id, None)
        if not gene_id:
            continue
        gene_to_tx_list_with_alu_peptides_dict.setdefault(gene_id, []).append(tx_id)

    for gene_id, checked_tx_list in gene_to_tx_list_with_alu_peptides_dict.items():
        gene_name = gene_id_to_name_dict.get(gene_id, gene_id)
        te_bed_set = set()
        for tx_id in checked_tx_list:
            bed_list_te = tx_matched_te_bed_dict.get(tx_id, None)
            te_bed_set |= set(bed_list_te)
            #if bed_list_te:
            #    break
        #else:
        #    continue #No Alu exons
        if len(te_bed_set) == 0:
            continue
        exon_dict = {tx_id:tx_exon_coord_dict[tx_id] for tx_id in checked_tx_list}
        checked_te_idx_dict = {}
        checked_te_list = []
        all_te_exon_list = []
        for arr in te_bed_set:
            checked_te_idx_dict[arr[3]] = checked_te_idx_dict.get(arr[3], 0) + 1
            te_idx = checked_te_idx_dict[arr[3]]
            te_id = '%s %d'%(arr[3], te_idx)
            checked_te_list.append(te_id)
            exon_dict.setdefault(te_id, []).append((int(arr[1]), int(arr[2])))
            all_te_exon_list.append((int(arr[1]), int(arr[2])))
        min_pos, max_pos = get_min_max_pos(exon_dict)
        final_pos_dict, corrected_exon_dict, corrected_intron_dict = adjust_pos(exon_dict, min_pos)
        checked_peptide_set = set()
        corrected_peptide_coord_count_dict = {}
        for tx_id in checked_tx_list:
            peptide_dict = tx_peptides_bed_dict.get(tx_id, {})
            for peptide_coord in peptide_dict.keys():
                if peptide_coord in checked_peptide_set:
                    continue
                peptide_start = int(peptide_coord[2])
                peptide_end = int(peptide_coord[3])
                if check_te_with_tx_exons((peptide_start, peptide_end), all_te_exon_list):
                    peptide_seq = peptide_dict[peptide_coord]
                    peptide_count = peptide_count_dict[peptide_seq]
                    peptide_start -= min_pos
                    peptide_end -= min_pos
                    corrected_peptide_coord_count_dict[(peptide_start, peptide_end)] = peptide_count
                    checked_peptide_set.add(peptide_coord)
        if len(corrected_peptide_coord_count_dict) > 0:
            total_peptide_count = sum(corrected_peptide_coord_count_dict.values())
            plot_structures(checked_tx_list, checked_te_list, corrected_exon_dict, final_pos_dict, 
                    peptide_count_dict=corrected_peptide_coord_count_dict,
                    output_file='%s/%s_PSM%d.png'%(output_folder, gene_name, total_peptide_count))



