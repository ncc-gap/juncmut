#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 09:23:51 2019
@author: genome
python juncmut_intersect.py <folder of query file> <file prefix>
"""

def juncmut_intersect(input_file, output_file, original_sj_file):
    import subprocess
    import os

    tmpfile1 = output_file + ".tmp1"
    tmpfile2 = output_file + ".tmp2"

    with open(tmpfile1, 'w') as hout:
        with open(input_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                print(F[0] + '\t' + str(int(F[1]) - 5) + '\t' + str(int(F[2]) + 5) + '\t' + '\t'.join(F), file = hout)

    hout = open(tmpfile2, 'w')
    subprocess.run(["bedtools", "intersect", "-a", tmpfile1, "-b", original_sj_file, "-c"], stdout = hout)
    hout.close()

    with open(output_file, 'w') as hout:
        with open(tmpfile2, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                print('\t'.join(F[3:]), file = hout)

             
    os.remove(tmpfile2)
    os.remove(tmpfile1)

if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--input_file", metavar = "input_file", default = None, type = str,
                            help = "input file") 
        
    parser.add_argument("--output_file", metavar = "output_file", default = "my_sample", type = str,
                            help = "output file") 
    
    parser.add_argument("--original_sj_file", metavar = "original_sj_file", default = "my_sample", type = str,
                            help = "original_sj_file") 
        
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    original_sj_file = args.original_sj_file
    
    juncmut_intersect(input_file, output_file, original_sj_file)

                
            
