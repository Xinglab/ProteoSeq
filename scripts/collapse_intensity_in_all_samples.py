import os,argparse,re
import numpy as np

def get_intensity_dict(input_pin):
    specid_intensity_dict = {}
    header_dict = {}
    intensity_column = -1
    log_base = 0
    for line in open(input_pin, 'r'):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('DefaultDirection'):
            continue
        if line.startswith('SpecId'):
            for i in range(len(arr)):
                if arr[i] == 'lnMS2IonCurrent':
                    intensity_column = i
                    log_base = 'exp'
                    break
                if arr[i] == 'log10_intensity':
                    intensity_column = i
                    log_base = 10
                    break
            continue
        specid = arr[0] #.split('_')[2]
        if intensity_column >= 0:
            intensity = float(arr[intensity_column])
            if log_base == 'exp':
                intensity = np.exp(intensity)
            elif log_base == 10:
                intensity = 10 ** intensity
        #specid_intensity_dict.setdefault(specid, set()).add(intensity) 
        specid_intensity_dict[specid] = intensity
    return specid_intensity_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collapse Percolator outputs to tabular table.')
    parser.add_argument('-i', '--input_file_list', help='List of input files with PSMs, space splitted. Required', nargs='+', type=str, required=True)
    parser.add_argument('-p', '--input_pin_list', help='List of input PIN files with MS2 spectra intensity, space splitted. Required', nargs='+', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=str)
    args = parser.parse_args()
    
    input_file_list = args.input_file_list
    input_pin_list = args.input_pin_list
    output_file = args.output_file

    if len(input_file_list) != len(input_pin_list):
        exit('Input psm files and input pin files do not have same length! Exit.')
    
    peptide_tx_dict = {}
    peptide_sample_intensity_dict = {}
    sample_set = set()
    for i in range(len(input_file_list)):
        input_pin = input_pin_list[i]
        specid_intensity_dict = get_intensity_dict(input_pin)
        input_file = input_file_list[i]
        sample = input_file.split('/')[-1].replace('.psms', '')
        sample_set.add(sample)
        for line in open(input_file, 'r'):
            arr = line.rstrip('\n').split('\t')
            if line.startswith('PSMId'):
                continue
            intensity = specid_intensity_dict[arr[0]]
            peptide_seq = arr[4].upper()
            #peptide_seq = re.sub('\[UNIMOD:\d+\]', '', peptide_seq.split('.')[1])
            peptide_seq = re.sub('\[.+?\]', '', peptide_seq)
            peptide_seq = re.findall('\.(.+)\.', peptide_seq)[0]
            peptide_sample_intensity_dict[peptide_seq][sample] = peptide_sample_intensity_dict.setdefault(peptide_seq, {}).get(sample, 0) + intensity
            tx_id_all = ';'.join(arr[5:])
            peptide_tx_dict.setdefault(peptide_seq, set()).add(tx_id_all)

    peptide_to_protein_dict = {}
    for peptide_seq,mapping_set in peptide_tx_dict.items():
        if len(mapping_set)>1:
            print('%s has multiple annotations:%s'%(peptide_seq, ' '.join(mapping_set)))
        for mapping in mapping_set:
            for m in mapping.split(';'):
                protein_idx = int(m.replace('Protein_', ''))
                peptide_to_protein_dict.setdefault(peptide_seq, set()).add(protein_idx)

    output = open(output_file, 'w')
    sample_list = sorted(sample_set)
    output.write('#peptide\tprotein_index\ttotal_intensity\t%s\n'%('\t'.join(sample_list)))
    for peptide_seq in sorted(peptide_sample_intensity_dict, key=lambda k:sum(peptide_sample_intensity_dict[k].values()), reverse=True):
        intensity_list = [peptide_sample_intensity_dict[peptide_seq].get(sample, 0) for sample in sample_list]
        total_intensity = sum(intensity_list)
        intensity_list = [('%.2f'%_).rstrip("0").rstrip(".") for _ in intensity_list]
        protein_idx_output = [str(_) if type(_)==int else _ for _ in sorted(peptide_to_protein_dict[peptide_seq])]
        output.write('%s\t%s\t%d\t%s\n'%(peptide_seq, ';'.join(protein_idx_output), total_intensity, '\t'.join(intensity_list)))        

    output.close()


