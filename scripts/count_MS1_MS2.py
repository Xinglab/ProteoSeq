import os,argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count raw spectrum in metadata files.')
    parser.add_argument('-i', '--input_folder', help='Input folder of metadata files. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table of count matrix. Required', type=str, required=True)
    args = parser.parse_args()

    input_folder = args.input_folder
    output_table = args.output_table

    f_list = os.listdir(input_folder)
    sample_ms1_count_dict = {}
    sample_ms2_count_dict = {}
    if len(sample_ms1_count_dict) == 0:
        for f in f_list:
            if not f.endswith('.txt'):
                continue
            sample = f.split('.')[0].replace('-metadata', '')
            for line in open('%s/%s'%(input_folder, f), 'r'):
                if line.find('Number of MS1 spectra') >= 0:
                    ms1 = re.findall('Number of MS1 spectra=\D*(\d+)', line)[0]
                    sample_ms1_count_dict[sample] = int(ms1)
                if line.find('Number of MS2 spectra') >= 0:
                    ms2 = re.findall('Number of MS2 spectra=\D*(\d+)', line)[0]
                    sample_ms2_count_dict[sample] = int(ms2)
    if len(sample_ms1_count_dict) == 0:
        flag = False
        for f in f_list:
            if not f.endswith('.txt'):
                continue
            sample = f.split('.')[0].replace('-metadata', '')
            for line in open('%s/%s'%(input_folder, f), 'r'):
                if line.find('Number of spectra per MS level:') >= 0:
                    flag = True
                elif line[0] != ' ':
                    flag = False
                elif flag:
                    if line.find('level 1') >= 0:
                        ms1 = re.findall('level 1: (\d+)', line.strip())[0]
                        sample_ms1_count_dict[sample] = int(ms1)
                    if line.find('level 2') >= 0:
                        ms2 = re.findall('level 2: (\d+)', line.strip())[0]
                        sample_ms2_count_dict[sample] = int(ms2)
    if len(sample_ms1_count_dict) == 0:
        exit('Wrong metadata format!')
    sample_list = sorted(sample_ms1_count_dict)
    output = open(output_table, 'w')
    output.write('#MS_Run\tMS1_Count\tMS2_Count\n')
    for sample in sample_list:
        output.write('%s\t%d\t%d\n'%(sample, sample_ms1_count_dict[sample], sample_ms2_count_dict[sample]))

    output.close()


