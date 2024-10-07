import os,argparse,re
import numpy as np
from scipy import stats,sparse,io

def multiply(a, b, block=10000):
    if (a.shape[1] != len(b)):
        print('Two arrays do not have matching shapes')
        return False
    nrow = a.shape[0]
    n = nrow//block
    temp_matrix = a[0:block, :].multiply(b)
    if n == 0:
        return temp_matrix.tocsr()
    result_list = []
    result_list.append(temp_matrix)
    i = 1
    while i < n:
        temp_matrix = a[i*block:(i+1)*block, :].multiply(b)
        result_list.append(temp_matrix)
        i += 1
    temp_matrix = a[i*block:nrow, :].multiply(b)
    result_list.append(temp_matrix)
    result = sparse.vstack(result_list)
    return result.tocsr()

def dot(a, b, block=10000):
    if a.shape[1] != len(b):
        print('Two arrays do not have matching shapes')
        return False
    nrow = a.shape[0]
    n = nrow//block
    temp_array = a[0:block, :].dot(b)
    if n == 0:
        return temp_array.reshape(nrow, 1)
    i = 1
    while i < n:
        temp_row = a[i*block:(i+1)*block, :].dot(b)
        temp_array = np.append(temp_array, temp_row)
        i += 1
    temp_row = a[i*block:nrow, :].dot(b)
    temp_array = np.append(temp_array, temp_row)
    return temp_array.reshape(nrow, 1)

def divide(a, b, block=10000):
    if (a.shape[0] != b.shape[0]):
        print('Two arrays do not have matching shapes')
        return False
    nrow = a.shape[0]
    n = nrow//block
    result_list = []
    temp_matrix = a[0:block, :].multiply(1 / b[0:block, :])
    if n == 0:
        return temp_matrix.tocsr()
    result_list.append(temp_matrix)
    i = 1
    while i < n:
        temp_matrix = a[i*block:(i+1)*block, :].multiply(1 / b[i*block:(i+1)*block, :])
        result_list.append(temp_matrix)
        i += 1
    temp_matrix = a[i*block:nrow, :].multiply(1 / b[i*block:nrow, :])
    result_list.append(temp_matrix)
    result = sparse.vstack(result_list)
    return result.tocsr()

def e_step(theta_t, block=10000):
    numerator = multiply(_Y, theta_t, block=block)
    denominator = dot(_Y, theta_t, block=block)
    z_t1 = divide(numerator, denominator, block=block) #numerator / denominator
    del numerator,denominator
    return z_t1

def m_step(z_t1):
    n_t1 = np.array(z_t1.sum(axis=0)).flatten()
    theta_t1 = n_t1/_N
    del n_t1
    return theta_t1

def print_dist(output_file, dist, idx):
    output_dist = open(output_file, 'a')
    output_dist.write('%s\t%f\n'%(idx, np.log(dist)))
    output_dist.close()

def print_theta(output_file, theta_t, idx):
    output_theta = open(output_file, 'a')
    output_theta.write('%d\t%s\n'%(idx, '\t'.join([str(_) for _ in theta_t])))
    output_theta.close()

def print_count(output_file, z_total):
    output_count = open(output_file, 'w')
    output_count.write('#%s\t%s\n'%('gene_name', 'peptide_count'))
    for _ in range(len(z_total)):
        gene_count = z_total[_]
        gene_name = all_mapped_gene_list[_]
        output_count.write('%s\t%.2f\n'%(gene_name, gene_count))
    output_count.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform EM algorithm to get gene level peptide abundance.')
    parser.add_argument('-i', '--index_file', help='Index file for protein sequences. Required', type=str, required=True)
    parser.add_argument('-p', '--psm_file', help='Input table of PSM count. Required', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output table of gene level expression. Required', type=str, required=True)
    parser.add_argument('-e', '--epsilon', help='Threshold for M step', type=float, default=0.0001)

    args = parser.parse_args()
    
    index_file = args.index_file
    psm_file = args.psm_file
    output_file = args.output_file

    epsilon = float(args.epsilon)

    protein_idx_all_gene_dict = {}
    header_dict = {}
    for line in open(index_file, 'r'):
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_dict[arr[i]] = i
            continue
        arr = line.rstrip('\n').split('\t')
        protein_idx = int(arr[header_dict['protein_idx']].replace('Protein_', ''))
        tx_id = arr[header_dict['seq_id']]
        if arr[header_dict['gene_name']] == '':
            continue
        gene_name_list = arr[header_dict['gene_name']].split(',')
        for gene_name in gene_name_list:
            if (gene_name != '') and (gene_name != 'NA'):
                protein_idx_all_gene_dict.setdefault(protein_idx, set()).add(gene_name)

    sample_list = []
    sample_peptide_count_dict = {}
    sample_peptide_mapped_idx_dict = {}
    header_rev_dict = {}
    for line in open(psm_file, 'r'):
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            arr = line[1:].rstrip('\n').split('\t')
            for i in range(len(arr)):
                header_rev_dict[i] = arr[i]
            sample_list = arr[2:]
            continue
        arr = line.rstrip('\n').split('\t')
        peptide_seq = arr[0]
        for i in range(2, len(arr)):
            sample = header_rev_dict[i]
            psms = float(arr[i])
            protein_idx_list = [int(_) for _ in arr[1].split(';')]
            if psms > 0:
                sample_peptide_count_dict.setdefault(sample, {})[peptide_seq] = psms
                sample_peptide_mapped_idx_dict.setdefault(sample, {})[peptide_seq] = protein_idx_list

    sample_gene_count_dict = {}
    for sample in sample_list:
        peptide_mapped_idx_dict = sample_peptide_mapped_idx_dict.get(sample, {})
        peptide_count_dict = sample_peptide_count_dict.get(sample, {})
        if len(peptide_mapped_idx_dict) == 0 or len(peptide_count_dict) == 0:
            continue
        peptide_mapped_gene_dict = {}
        for peptide_seq, protein_idx_list in peptide_mapped_idx_dict.items():
            for protein_idx in protein_idx_list:
                gene_set = protein_idx_all_gene_dict.get(protein_idx, None)
                if not gene_set:
                    continue
                peptide_mapped_gene_dict[peptide_seq] = peptide_mapped_gene_dict.setdefault(peptide_seq, set()) | gene_set

        peptide_mapped_gene_dict_for_training = peptide_mapped_gene_dict.copy()
            
        all_mapped_gene_idx_dict = {}
        idx = 0
        for peptide_seq, mapped_gene_set in peptide_mapped_gene_dict_for_training.items():
            for gene_name in mapped_gene_set:
                if (gene_name=='') or (gene_name=='NA'):
                    continue
                gene_idx = all_mapped_gene_idx_dict.get(gene_name, -1)
                if gene_idx < 0:
                    all_mapped_gene_idx_dict[gene_name] = idx
                    idx += 1

        all_mapped_gene_list = sorted(all_mapped_gene_idx_dict, key=lambda k:all_mapped_gene_idx_dict[k])

        peptide_idx_dict = {}
        peptide_idx_count_dict = {}
        idx = 0
        for peptide_seq in sorted(peptide_mapped_gene_dict_for_training):
            peptide_idx_dict[peptide_seq] = idx
            peptide_idx_count_dict[idx] = peptide_count_dict[peptide_seq]
            idx += 1
        peptide_mapped_gene_idx_dict = {}
        for peptide_seq, gene_set in peptide_mapped_gene_dict_for_training.items():
            gene_idx_set = set(all_mapped_gene_idx_dict[_] for _ in gene_set)
            peptide_idx = peptide_idx_dict[peptide_seq]
            peptide_mapped_gene_idx_dict[peptide_idx] = gene_idx_set


        row = []
        column = []
        psm_count = []
        for i in sorted(peptide_mapped_gene_idx_dict):
            psm_count.append(peptide_idx_count_dict[i])
            j_set = peptide_mapped_gene_idx_dict[i]
            row += [i]*len(j_set)
            column += list(j_set)

        elements = [1]*len(row)
        _Y = sparse.csr_array((elements, (row, column)))
        _N, _K = _Y.shape
        #print("Size of array:%dX%d"%(_N, _K))
        
        #io.mmwrite('%s/%s_indicator_matrix'%(output_folder, output_prefix), _Y)
        #io.mmwrite('%s_indicator_matrix'%(output_prefix), _Y)

        del row, column, elements, peptide_mapped_gene_idx_dict

        theta_t = np.array([1/_K]*_K)

        i = 0
        #theta_file = '%s/%s_gene_count_estimated_by_em_theta.txt'%(output_folder, output_prefix)
        #theta_file = '%s_gene_count_estimated_by_em_theta.txt'%(output_prefix)
        #output_theta = open(theta_file, 'w')
        #output_theta.write('#%s\t%s\n'%('iteration', '\t'.join(all_mapped_gene_list)))
        #theta_line = '%d\t%s\n'%(i, '\t'.join([str(_) for _ in theta_t])) #initial probability, 0
        #output_theta.write(theta_line)
        #output_theta.close()

        #dist_file = '%s/%s_gene_count_estimated_by_em_dist.txt'%(output_folder, output_prefix)
        #dist_file = '%s_gene_count_estimated_by_em_dist.txt'%(output_prefix)
        #output_dist = open(dist_file, 'w')
        #output_dist.write('#%s\t%s\n'%('iteration', 'log_dist'))
        #output_dist.close()

        max_iter = 10000
        while i < max_iter:
            i += 1
            z_t1 = e_step(theta_t, block=50000)
            theta_t1 = m_step(z_t1)
            dist = np.sum(abs(theta_t1 - theta_t))
            #print_dist(dist_file, dist, i) #output to existing file
            theta_t = theta_t1
            del theta_t1
            #print_theta(theta_file, theta_t, i)
            #z_total = np.array(z_t1.sum(axis=0)).flatten()
            z_t1_wpsm = multiply(z_t1.transpose(), psm_count).transpose()
            z_total = np.array(z_t1_wpsm.sum(axis=0)).flatten()
            if dist < epsilon:
                for _ in range(len(z_total)):
                    gene_count = z_total[_]
                    gene_name = all_mapped_gene_list[_]
                    sample_gene_count_dict.setdefault(gene_name, {})[sample] = gene_count
                break
            del z_t1
        
        if i >= max_iter:
            print('Does not converge after %d times iteration in sample: %s'%(max_iter, sample))
            for _ in range(len(z_total)):
                gene_count = z_total[_]
                gene_name = all_mapped_gene_list[_]
                sample_gene_count_dict.setdefault(gene_name, {})[sample] = gene_count
    
    output_count = open(output_file, 'w')
    output_count.write('#%s\t%s\n'%('gene_name', '\t'.join(sample_list)))
    for gene_name in sorted(sample_gene_count_dict):
        gene_count_dict = sample_gene_count_dict[gene_name]
        gene_count_list = [str(gene_count_dict.get(_,0)) for _ in sample_list]
        output_count.write('%s\t%s\n'%(gene_name, '\t'.join(gene_count_list)))
    output_count.close()




