"""
Script Name: plot_peptides_to_tx.py
Description: Plot peptides to the transcripts structures
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os,re,sys,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
from matplotlib import ticker
from matplotlib.patches import Rectangle

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


def get_min_pos(exon_coord_list, pos_neg='+'):
    if pos_neg == '+':
        min_pos = 1000000000
        for exon_coord in exon_coord_list:
            min_pos = min(min_pos, min(exon_coord))
    elif pos_neg == '-':  #negative values
        min_pos = -1000000000
        for exon_coord in exon_coord_list:
            min_pos = max(min_pos, max(exon_coord))
    return min_pos


def get_max_pos(exon_coord_list, pos_neg='+'):
    if pos_neg == '+':
        max_pos = 0
        for exon_coord in exon_coord_list:
            max_pos = max(max_pos, max(exon_coord))
    elif pos_neg == '-': #negative values
        max_pos = 0
        for exon_coord in exon_coord_list:
            max_pos = min(max_pos, min(exon_coord))
    return max_pos


def get_tss_tes_pos(exon_dict, strand):
    concat_coord_list = []
    for tx_id, exon_coord_list in exon_dict.items():
        concat_coord_list += exon_coord_list
    if strand == '+':
        tss_pos = get_min_pos(concat_coord_list, pos_neg='+')
        tes_pos = get_max_pos(concat_coord_list, pos_neg='+')
    elif strand == '-':
        tss_pos = get_max_pos(concat_coord_list, pos_neg='+')
        tes_pos = get_min_pos(concat_coord_list, pos_neg='+')
    return tss_pos, tes_pos


def adjust_pos(exon_coord, tss_pos):
    return (exon_coord[0]-tss_pos, exon_coord[1]-tss_pos)


def adjust_pos_all_exon(exon_coord_list, tss_pos):
    if len(exon_coord_list) == 0:
        return []
    adjusted_exon_coord_list = []
    for exon_coord in exon_coord_list:
        adjusted_exon_coord = adjust_pos(exon_coord, tss_pos)
        adjusted_exon_coord_list.append(adjusted_exon_coord)
    return adjusted_exon_coord_list


def adjust_pos_all_tx(exon_dict, tss_pos):
    adjusted_exon_dict = {}
    for tx_id, exon_coord_list in exon_dict.items():
        adjusted_exon_coord_list = adjust_pos_all_exon(exon_coord_list, tss_pos)
        adjusted_exon_dict[tx_id] = adjusted_exon_coord_list
    return adjusted_exon_dict


def select_color(tx_id):
    if tx_id.startswith('ENST'):
        return 'blue'
    return 'red'


def select_hatch(tx_id):
    if tx_id in peptide_tx_list:
        return '/'
    return None


def draw_intron(ax, y_pos, exon_coord_list, strand, color='k', linewidth=1, linestyle='-'):
    if len(exon_coord_list) < 2:
        return False
    intron_left = exon_coord_list[0][1]
    intron_right = exon_coord_list[-1][0]
    if strand == '+':
        ax.hlines(y_pos, intron_left, intron_right, color=color, lw=linewidth, zorder=1, ls=linestyle)
    else:
        ax.hlines(y_pos, -intron_right, -intron_left, color=color, lw=linewidth, zorder=1, ls=linestyle)


def draw_single_exon(ax, y_pos, exon_coord, max_pos, strand, color='#D9D9D9', exon_height=.6):
    (exon_left, exon_right) = exon_coord
    exon_length = exon_right - exon_left + 1 #????
    if exon_length < max_pos*0.001:
        exon_length = max_pos*0.001
    if strand == '+':
        ax.add_patch(Rectangle((exon_left, y_pos-.5*exon_height), exon_length, exon_height, fill=True, fc=color, lw=exon_border, edgecolor='k'))
    else:
        ax.add_patch(Rectangle((-exon_right, y_pos-.5*exon_height), exon_length, exon_height, fill=True, fc=color, lw=exon_border, edgecolor='k'))


def draw_all_exon(ax, y_pos, exon_coord_list, max_pos, strand, cds_coord_list=[], 
                  exon_color='#D9D9D9', exon_height=.6, plot_cds=False):
    if len(exon_coord_list) == 0:
        return False
    exon_coord_list = sorted(exon_coord_list)
    for exon_coord in exon_coord_list:
        if plot_cds and (exon_coord in cds_coord_list):
            continue
        draw_single_exon(ax, y_pos, exon_coord, max_pos, strand, color=exon_color, exon_height=exon_height)
    if plot_cds:
        cds_coord_list = sorted(cds_coord_list)
        for exon_coord in cds_coord_list:
            draw_single_exon(ax, y_pos, exon_coord, max_pos, strand, color=exon_color, exon_height=exon_height*1.5)


def plot_structures(ax, tx_list, exon_coord_dict, tx_id_label_list=None, cds_coord_dict={}, plot_cds_in_tx=True, peptide_coord_dict={}, exon_height=.6, intron_line_weight=1, peptide_junc_weight=2):
    num_trans = len(tx_list)
    peptide_list = sorted(peptide_coord_dict, key=lambda k:peptide_coord_dict[k])
    num_peptides = len(peptide_list)
    ymin = -num_peptides + .5
    ymax = num_trans+.5
    ax.set_ylim(ymin, ymax)
    max_final_pos = 0
    for tx_idx in range(num_trans):
        tx_id = tx_list[tx_idx]
        strand = tx_strand_dict[tx_id]
        exon_coord_list = exon_coord_dict[tx_id]                           
        max_pos = get_max_pos(exon_coord_list, pos_neg=strand) #adjusted, negative valaues for antisense
        max_final_pos = get_max_pos([(max_final_pos, max_pos)], pos_neg=strand)
    if strand == '-':
        max_final_pos = -max_final_pos
        peptide_list = peptide_list[::-1]
    xmin, xmax = -max_final_pos*.01, max_final_pos*1.01
    for tx_idx in range(num_trans):
        y_pos = num_trans - tx_idx
        tx_id = tx_list[tx_idx]
        if tx_id.startswith('ENST'):
            exon_color = '#4472C4'
        else:
            exon_color = '#FF5333'
        strand = tx_strand_dict[tx_id]
        exon_coord_list = exon_coord_dict[tx_id]
        cds_coord_list = cds_coord_dict.get(tx_id, [])
        draw_intron(ax, y_pos, exon_coord_list, strand, 'k', intron_line_weight, '-')
        draw_all_exon(ax, y_pos, exon_coord_list, max_final_pos, strand, cds_coord_list=cds_coord_list, plot_cds=plot_cds_in_tx, exon_color=exon_color, exon_height=exon_height)
    for peptide_idx in range(num_peptides):
        #y_pos = peptide_idx - num_peptides + 1
        y_pos = - peptide_idx
        peptide_seq = peptide_list[peptide_idx]
        peptide_coord_list = peptide_coord_dict[peptide_seq]
        draw_intron(ax, y_pos, peptide_coord_list, strand, '#008000', peptide_junc_weight, '--')
        draw_all_exon(ax, y_pos, peptide_coord_list, max_final_pos, strand, exon_color='#008000', exon_height=exon_height)
    ax.set_xlim(xmin, xmax)
    for tx_idx in range(num_trans):
        y_pos = num_trans - tx_idx
        tx_id = tx_list[tx_idx]
        strand = tx_strand_dict[tx_id]
        if tx_id_label_list and len(tx_id_label_list) >= num_trans:
            tx_id_label = tx_id_label_list[tx_idx]
        else:
            tx_id_label = tx_id
        if strand == '+':
            #ax.text(0, y_pos+exon_height/2, tx_id_label, va='bottom', ha='left', fontsize=font_size+1)
            ax.text(xmin, y_pos, tx_id_label, va='center', ha='right', fontsize=font_size)
        else:
            #ax.text(-max_final_pos, y_pos+exon_height/2, tx_id_label, va='bottom', ha='left', fontsize=font_size+1)
            ax.text(xmin, y_pos, tx_id_label, va='center', ha='right', fontsize=font_size)
    for peptide_idx in range(num_peptides):
        #y_pos = peptide_idx - num_peptides + 1
        y_pos = - peptide_idx
        peptide_seq = peptide_list[peptide_idx]
        if strand == '+':
            ax.text(xmin, y_pos, 'Peptide_%s'%peptide_seq, va='center', ha='right', fontsize=font_size)
        else:
            ax.text(xmin, y_pos, 'Peptide_%s'%peptide_seq, va='center', ha='right', fontsize=font_size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks([], minor=False)
    ax.set_xticks([0, max_final_pos])
    tss_label = '{:,d}'.format(tss_pos) # '%.1f'%(tss_pos/1000000)
    tes_label = '{:,d}'.format(tss_pos+max_final_pos) #'%.1f'%((tss_pos+max_final_pos)/1000000)
    ax.set_xticklabels([tss_label, tes_label])
    chrom = '%s'%tx_chrom_dict[tx_id_list[0]] # (Mb)
    ax.text(.5, -.02, chrom, va='top', ha='center', fontsize=font_size, transform = ax.transAxes)
    #ax.set_xlabel(chrom, fontsize=font_size+2, labelpad=-10)
    #ax.axis('off')


def get_tx_id_label_short(tx_id):
    canonical = tx_canonical_dict.get(tx_id, None)
    if canonical:
        return '%s (Canonical)'%tx_id
    tx_id_label_list = []
    orf_type = tx_orf_type_dict.get(tx_id, None)
    if orf_type == 'GENCODE':
        tx_id_label_list.append('ORF:GENCODE annotated')
    elif orf_type == 'novel':
        tx_id_label_list.append('ORF:predicted')
    #else:
    #    tx_id_label_list.append('Not in protein db')
    if len(tx_id_label_list) >= 1:
        return '%s (%s)'%(tx_id, '; '.join(tx_id_label_list))
    return '%s'%(tx_id)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot peptides to the transcripts structures.')
    parser.add_argument('-i', '--input_gtf', help='Input GTF file of peptide mappings. Required', type=str, required=True)
    parser.add_argument('-r', '--ref_gtf', help='Input GTF file of protein references. Required', type=str, required=True)
    parser.add_argument('-o', '--output_folder', help='Output folder to plot figures', type=str, required=True)
    parser.add_argument('--exon_height', help='Exon height, default=.3', default=0.3, type=float)
    parser.add_argument('--exon_border', help='Exon border, default=0', default=0, type=float)
    parser.add_argument('--intron_line_weight', help='Intron line weight, default=.5', default=0.5, type=float)
    parser.add_argument('--peptide_junc_weight', help='Peptide junction dashed line weight, default=1', default=1, type=float)
    parser.add_argument('--font_size', help='Font size, default=10', default=10, type=int)
    args = parser.parse_args()

    input_gtf = args.input_gtf
    ref_gtf = args.ref_gtf
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    exon_height = float(args.exon_height)
    exon_border = float(args.exon_border)
    intron_line_weight = float(args.intron_line_weight)
    peptide_junc_weight = float(args.peptide_junc_weight)
    font_size = float(args.font_size)

    gene_id_to_name_dict = {}
    detected_gene_id_to_tx_dict = {}
    detected_tx_to_gene_id_dict = {}
    tx_coord_dict = {}
    tx_exon_coord_dict = {}
    tx_chrom_dict = {}
    tx_strand_dict = {}
    tx_orf_coord_dict = {}
    tx_orf_type_dict = {}
    tx_canonical_dict = {}
    for line in open(ref_gtf, 'r'):
        if line[0] == '#':
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'transcript':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            tx_coord = (int(arr[3]), int(arr[4]))
            tx_coord_dict[tx_id] = tx_coord
            gene_id = d.get('gene_id', [''])[0]
            gene_name = d.get('gene_name', [''])[0]
            if gene_id != '':
                detected_tx_to_gene_id_dict[tx_id] = gene_id
                detected_gene_id_to_tx_dict.setdefault(gene_id, []).append(tx_id)
                if gene_name != '':
                    gene_id_to_name_dict[gene_id] = gene_name
            tx_chrom_dict[tx_id] = arr[0]
            tx_strand_dict[tx_id] = arr[6]
        elif arr[2] == 'exon':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            exon_coord = (int(arr[3]), int(arr[4]))
            tx_exon_coord_dict.setdefault(tx_id, []).append(exon_coord)
            orf_type = d.get('ORF', [''])[0]
            if orf_type != '':
                tx_orf_type_dict[tx_id] = orf_type
            if 'Ensembl_canonical' in d.get('tag', []):
                tx_canonical_dict[tx_id] = 1
        elif arr[2] == 'CDS':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            cds_coord = (int(arr[3]), int(arr[4]))
            tx_orf_coord_dict.setdefault(tx_id, []).append(cds_coord)

    gene_peptide_seq_coord_dict = {}
    for line in open(input_gtf, 'r'):
        if line[0] == '#':
            continue
        arr = line.rstrip('\n').split('\t')
        if arr[2] == 'peptide':
            d = parse_attributes(arr[8])
            tx_id = d.get('transcript_id', [''])[0]
            if tx_id == '':
                continue
            gene_id = detected_tx_to_gene_id_dict.get(tx_id, '')
            if gene_id == '' or gene_id == 'NA' or gene_id == 'None':
                continue
            peptide_id = d['peptide_id'][0]
            peptide_seq = d['peptide_seq'][0]
            peptide_coord = (int(arr[3]), int(arr[4]))
            gene_peptide_seq_coord_dict.setdefault(gene_id, {}).setdefault(tx_id, {}).setdefault(peptide_seq, []).append(peptide_coord)

    #for peptide_tx_id in tx_peptide_seq_coord_dict.keys():
    for gene_id in gene_peptide_seq_coord_dict.keys():
        tx_peptide_seq_coord_dict = gene_peptide_seq_coord_dict[gene_id]
        peptide_tx_list = sorted(tx_peptide_seq_coord_dict.keys())
        peptide_tx_id = peptide_tx_list[0]
        peptide_seq_coord_dict = tx_peptide_seq_coord_dict[peptide_tx_id]
        gene_name = gene_id_to_name_dict.get(gene_id, gene_id)
        tx_id_list = detected_gene_id_to_tx_dict[gene_id]
        tx_id_label_list = []
        for tx_id in tx_id_list:
            tx_id_label = get_tx_id_label_short(tx_id)
            tx_id_label_list.append(tx_id_label)
        strand = tx_strand_dict[tx_id_list[0]]
        exon_dict = {tx_id:tx_exon_coord_dict[tx_id] for tx_id in tx_id_list}
        tss_pos, tes_pos = get_tss_tes_pos(exon_dict, strand)
        adjusted_exon_dict = adjust_pos_all_tx(exon_dict, tss_pos)
        cds_dict = {tx_id:tx_orf_coord_dict.get(tx_id, set([])) for tx_id in tx_id_list}
        adjusted_cds_dict = adjust_pos_all_tx(cds_dict, tss_pos)
        adjusted_peptide_dict = adjust_pos_all_tx(peptide_seq_coord_dict, tss_pos)
        fig_height = 1 + .5*exon_height*(len(tx_id_list)+len(peptide_seq_coord_dict))
        fig, ax = plt.subplots(1, 1, figsize=(9, fig_height), sharey=True)
        plot_structures(ax, tx_id_list, adjusted_exon_dict, 
                tx_id_label_list=tx_id_label_list, cds_coord_dict=adjusted_cds_dict, plot_cds_in_tx=True,
                peptide_coord_dict=adjusted_peptide_dict, exon_height=exon_height,
                intron_line_weight=intron_line_weight, peptide_junc_weight=peptide_junc_weight)
        fig.tight_layout()
        fig.savefig('%s/%s_structure.pdf'%(output_folder, gene_name))
        plt.close()
        plt.clf()
                    

