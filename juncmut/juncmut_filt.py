#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Dec 2020

@author: naoko iida


Filtering
. select True.
. Overlap SJ <= 5
. remove IG gene
. gnomad <= 0.01
. ==> output 
"""

import os
import csv

def juncmut_filt(input_file, output_file):
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        header = hin.readline()
        hout.write(header)
        col = header.rstrip('\n').split('\t')
        for line in hin:
            F = line.rstrip('\n').split('\t')
            
            mut_key = F[col.index('Mut_key')].split(',')
            mut_key_chr = mut_key[0]
            mut_key_start = int(mut_key[1])

            sj_overlap_count = int(F[col.index('SJ_overlap_count')])
            gnomad_af = float(F[col.index('gnomAD_AF')])

            ## remove IG region
            #chr14:105586437..106879844
            if mut_key_chr == 'chr14' and mut_key_start >= 105586437 and mut_key_start <= 106879844:
                continue
            #chr2:88857361..90235368
            elif mut_key_chr == 'chr2' and mut_key_start >= 88857361 and mut_key_start <= 90235368:
                continue
            #chr22:22026076..22922913
            elif mut_key_chr == 'chr22' and mut_key_start >= 22026076 and mut_key_start <= 22922913:
                continue

            if sj_overlap_count <= 5 and gnomad_af <= 0.01:
                hout.write(line)

if __name__ == "__main__":
    import sys

    input_file = sys.argv[0]
    output_file = sys.argv[1]

    juncmut_filt(input_file, output_file)
