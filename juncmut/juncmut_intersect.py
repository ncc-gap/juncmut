#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 09:23:51 2019
@author: genome
python juncmut_intersect.py <folder of query file> <file prefix>
"""

def juncmut_intersect(pr, folder):
    import subprocess
    import pathlib
    import os
    
    fname = pr + '.SJ.fil.annot.assadjunifreqT.pmut.txt' 
    f = pathlib.Path(fname)
    fstem = f.stem
    
    input_file = './data/%s/alterativeSJ_mutprediction/%s' %(folder,fname)  #file = './alterativeSJ_cmut/%s/%s.SJ.fil.annot.assadjunifreqT.AI.5mut.txt'
    output_file = './data/%s/alterativeSJ_mutprediction/%s.SJinSJ.txt' %(folder,fstem)
    
    tmp1 = './data/%s/alterativeSJ_mutprediction/%s_tmp.bed' %(folder, pr)
    tmp2 = './data/%s/alterativeSJ_mutprediction/%s_tmp_intersect.bed' %(folder, pr)
    
    with open(input_file, 'r') as fin:
        with open(output_file, 'w') as fout:
            for line in fin:
                l = line.rstrip('\n')
                F = line.rstrip('\n').split('\t')
                sample = F[5]
                #print(sample)
                if(F[0]=="chr"):
                    rec = l + "\t" + "SJinSJcount\n"
                    fout.write(rec)
                else:
                    with open(tmp1, 'w') as tout:
                        position = F[0] + "\t" + str(int(F[1])-5) + "\t" + str(int(F[2])+5) + "\n"
                        tout.write(position)
                    
                    with open(tmp2, 'w') as iout:
                        path = "bedtools"
                        j_file = "./data/%s/%s.SJ.out.tab" %(folder, sample)
                        #fraction = str(1.0)
                        subprocess.run([path, "intersect", "-a", tmp1, "-b", j_file], stdout = iout) #"-loj",,  "-F", fraction
                    
                    num_lines = sum(1 for line in open(tmp2))
                    rec = l + "\t" + str(num_lines) + "\n"
                    fout.write(rec)
                    
    os.remove(tmp2)
    os.remove(tmp1)

if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input", metavar = "prefix", default = None, type = str,
                            help = "prefix of input file") 
        
    parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                            help = "Path to input file") 
        
    args = parser.parse_args()
    
    pr = args.input
    folder = args.folder #folder of query files
    
    juncmut_intersect(pr, folder)

                
            
