#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 09:23:51 2019
@author: genome
python juncmut_intersect.py <folder of query file> <file prefix>
"""

def juncmut_intersect(input_file, output_file, original_sj_file):
    import subprocess
    import pathlib
    import os

    # sample = os.path.basename(input_SJ).replace(".SJ.out.tab", '')    
    # fname = pr + '.SJ.fil.annot.assadjunifreqT.pmut.txt' 
    # f = pathlib.Path(fname)
    # fstem = f.stem
    
   #  input_file = './data/%s/alterativeSJ_mutprediction/%s' %(folder,fname)  #file = './alterativeSJ_cmut/%s/%s.SJ.fil.annot.assadjunifreqT.AI.5mut.txt'
    # output_file = './data/%s/alterativeSJ_mutprediction/%s.SJinSJ.txt' %(folder,fstem)
    
    # tmp1 = './data/%s/alterativeSJ_mutprediction/%s_tmp.bed' %(folder, pr)
    # tmp2 = './data/%s/alterativeSJ_mutprediction/%s_tmp_intersect.bed' %(folder, pr)
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

    """
    with open(input_file, 'r') as fin:
        with open(output_file, 'w') as fout:
            for line in fin:
                # l = line.rstrip('\n')
                F = line.rstrip('\n').split('\t')
                # sample = F[5]
                #print(sample)
                if(F[0]=="chr"):
                    rec = l + "\t" + "SJinSJcount\n"
                    fout.write(rec)
                else:
                    with open(tmpfile1, 'w') as tout:
                        position = F[0] + "\t" + str(int(F[1])-5) + "\t" + str(int(F[2])+5) + "\n"
                        tout.write(position)
                    
                    with open(tmpfile2, 'w') as iout:
                        path = "bedtools"
                        #$ j_file = "./data/%s/%s.SJ.out.tab" %(folder, sample)
                        #fraction = str(1.0)
                        subprocess.run([path, "intersect", "-a", tmpfile1, "-b", original_sj_file], stdout = iout) #"-loj",,  "-F", fraction
                    
                    num_lines = sum(1 for line in open(tmpfile2))
                    rec = l + "\t" + str(num_lines) + "\n"
                    fout.write(rec)
    """
             
    os.remove(tmpfile2)
    os.remove(tmpfile1)

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

                
            
