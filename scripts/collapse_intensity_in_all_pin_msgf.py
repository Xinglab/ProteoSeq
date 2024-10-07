import os,argparse,re
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collapse PIN outputs to tabular table.')
    parser.add_argument('-i', '--input_pin_list', help='List of input PIN files with MS2 spectra intensity, space splitted. Required', nargs='+', type=str, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=str)
    args = parser.parse_args()
    
    input_pin_list = args.input_pin_list
    output_file = args.output_file

    peptide_tx_dict = {}
    peptide_sample_intensity_dict = {}
    sample_set = set()
    for input_pin in input_pin_list:
        sample = input_pin.split('/')[-1].replace('.sig.pin', '')
        sample = sample.replace('.pin', '')
        sample_set.add(sample)
        header_dict = {}
        intensity_column = -1
        log_base = 0
        for line in open(input_pin, 'r'):
            arr = line.rstrip('\n').split('\t')
            if line.startswith('DefaultDirection'):
                continue
            if len(header_dict) == 0:
                for i in range(len(arr)):
                    header_dict[arr[i]] = i
                    if arr[i] == 'Ms1.Area':
                        intensity_column = i
                        log_base = 0
                    if arr[i] == 'lnMS2IonCurrent':
                        intensity_column = i
                        log_base = 'exp'
                    if arr[i] == 'log10_intensity':
                        intensity_column = i
                        log_base = 10
                continue
            if intensity_column >= 0:
                intensity = float(arr[intensity_column])
                if log_base == 'exp':
                    intensity = np.exp(intensity)
                elif log_base == 10:
                    intensity = 10 ** intensity
            else:
                exit('No intensity column is found in the pin file %s'%input_pin)
            peptide_seq = arr[header_dict['Peptide']].upper()
            peptide_seq = re.sub('\[.+?\]', '', peptide_seq)
            peptide_seq = re.findall('\.(.+)\.', peptide_seq)[0]
            tx_id_all = ';'.join(sorted(arr[header_dict['Proteins']:]))
            if tx_id_all.find('Random_') >= 0:
                continue
            peptide_sample_intensity_dict[peptide_seq][sample] = peptide_sample_intensity_dict.setdefault(peptide_seq, {}).get(sample, 0) + intensity 
            peptide_tx_dict.setdefault(peptide_seq, set()).add(tx_id_all)

    peptide_to_protein_dict = {}
    for peptide_seq, mapping_set in peptide_tx_dict.items():
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


