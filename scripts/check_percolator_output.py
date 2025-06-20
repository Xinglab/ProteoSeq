"""
Script Name: check_percolator_output.py
Description: If no significant PSM is returned by Percolator, check the log file to generate empty file if necessary.
Author: Lingyu Guan
Affiliation: Children's Hospital of Philadelphia (CHOP), Xing Lab
Email: guanl@chop.com
Date: 2025-06-19
"""

import os, argparse, re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='If no significant PSM is returned by Percolator, check the log file to generate empty file if necessary.')
    parser.add_argument('-i', '--percolator_output', help='Input output of Percolator. Required', type=str, required=True)
    parser.add_argument('-l', '--percolator_log', help='Input log of Percolator. Required', type=str, required=True)
    parser.add_argument('-o', '--output', help='Output adjusted output of Percolator. Required', type=str, required=True)
    args = parser.parse_args()
    
    percolator_output = args.percolator_output
    percolator_log = args.percolator_log
    output_file = args.output
    
    if os.path.exists(percolator_output):
        output = open(output_file, 'w')
        for line in open(percolator_output, 'r'):
            output.write(line)
        output.close()
    else:
        if not os.path.exists(percolator_log):
            print('Percolator failed! No output generated.')
        else:
            for line in open(percolator_log,'r'):
                if line.find('No targets found') >= 0:
                 output = open(output_file, 'w')
                 output.close()
                 print('Percolator is successfully run,  whereas no significant targets found. Touch a empty output for downstream process.')


