"""
Script Name: plot_peptides_spectra.py
Description: Plot MS/MS spectra matched to given peptide sequences.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os, sys, argparse
import yaml
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot MS/MS spectra matched to given peptide sequences.')
    parser.add_argument('-i', '--input_file', help='Input file of the peptides table. Required', type=str, required=True)
    parser.add_argument('-rd', '--result_dir', help='Input directory of peptides table of each dataset. Required', type=str, required=True)
    parser.add_argument('-wd', '--work_dir', help='Input directory of mzid files each dataset. Required', type=str, required=True)
    parser.add_argument('-c', '--config_file', help='YAML config file including the path to mzML directories. Required', type=str, required=True)
    parser.add_argument('-o', '--output_dir', help='Output directory. Required', type=str, required=True)
    args = parser.parse_args()

    input_file = args.input_file
    result_dir = args.result_dir
    work_dir = args.work_dir
    config_file = args.config_file
    output_dir = args.output_dir

    header_dict = {}
    dataset_peptides_dict = {}
    for line in open(input_file, 'r'):
        arr = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            for i in range(5, len(arr)):
                header_dict[i] = arr[i]
            continue
        for i in range(5, len(arr)):
            if int(arr[i]) > 0:
                dataset = header_dict[i]
                dataset_peptides_dict.setdefault(dataset, set()).add(arr[0])

    dataset_selected_sample_peptides_dict = {}
    for dataset, peptide_set in dataset_peptides_dict.items():
        dataset_input = f'{result_dir}/{dataset}.sig_psms_count.txt'
        for line in open(dataset_input, 'r'):
            arr = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                for i in range(3, len(arr)):
                    header_dict[i] = arr[i]
                continue
            if arr[0] not in peptide_set:
                continue
            for i in range(3, len(arr)):
                if int(arr[i]) > 0:
                    sample = header_dict[i]
                    dataset_selected_sample_peptides_dict.setdefault(dataset, {}).setdefault(sample, set()).add(arr[0])

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    mzml_dirs = {}
    for dataset in config.get('ms_datasets', []):
        mzml_dirs[dataset] = config['ms_datasets'][dataset].get('mzml_dir')

    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)

    importr('MSnbase')
    for dataset, sample_peptides_dict in dataset_selected_sample_peptides_dict.items():
        dataset_mzml_dir = mzml_dirs.get(dataset, '')
        for sample, peptide_set in sample_peptides_dict.items():
            input_mzid = f'{work_dir}/{dataset}/{sample}.mzid'
            input_mzml = f'{dataset_mzml_dir}/{sample}.mzML'
            if not os.path.exists(input_mzid):
                print(f'Error: {input_mzid} not exists!')
                continue
            if not os.path.exists(input_mzml):
                print(f'Error: {input_mzml} not exists!')
                continue
            ro.globalenv['pep_seq_list'] = StrVector(list(peptide_set))
            print(f'Processing dataset: {dataset}, sample: {sample}, peptides: {len(peptide_set)}')
            r(f"""
                iddata <- readMzIdData("{input_mzid}")
                rawdata <- readMSData("{input_mzml}", msLevel = 2, mode = "onDisk", verbose = FALSE)
                rawdata_wid <- addIdentificationData(rawdata, id = iddata, decoy = "isDecoy", rank = "rank", acc = "DatabaseAccess", desc = "DatabaseDescription", icol = "acquisitionNum", fcol = "acquisitionNum", pepseq = "sequence", accession = NULL, verbose = TRUE)
                for (pep_seq in pep_seq_list) {{
                for (i in which(grepl(pep_seq, fData(rawdata_wid)$sequence))) {{
                out_prefix <- paste(pep_seq, "{dataset}", "{sample}", i, sep = "_")
                out_pdf <- paste(out_prefix, "pdf", sep = ".")
                out_pdf_path <- paste("{output_dir}", out_pdf, sep = "/")
                pdf(out_pdf_path, width=7, height=5)
                plot(rawdata_wid[[i]], pep_seq, main = pep_seq)
                dev.off()
                }}
                }}
                rm(iddata, rawdata, rawdata_wid)
                gc()
            """)
    print(f'Plots saved to {output_dir}')
    print('Done!')