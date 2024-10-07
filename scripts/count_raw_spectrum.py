import os,argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count raw spectra in mzML files.')
    parser.add_argument('-i', '--input_folder', help='Input folder of mzML files. Required', type=str, required=True)
    parser.add_argument('-o', '--output_table', help='Output table of count matrix. Required', type=str, required=True)
    args = parser.parse_args()

    input_folder = args.input_folder
    output_table = args.output_table

    f_list = os.listdir(input_folder)
    sample_count_dict = {}
    for f in f_list:
        if (not f.endswith('.mzML')) and (not f.endswith('.mzml')):
            continue
        sample = f.split('.')[0]
        for line in open('%s/%s'%(input_folder, f), 'r'):
            res = re.findall('<spectrumList count=\"(\d+)\"', line)
            if len(res) > 0:
                count = int(res[0])
                sample_count_dict[sample] = count
                break
                
    sample_list = sorted(sample_count_dict)
    output = open(output_table, 'w')
    output.write('#MS_Run\tSpectrum_Count\n')
    for sample in sample_list:
        output.write('%s\t%d\n'%(sample, sample_count_dict[sample]))

    output.close()
