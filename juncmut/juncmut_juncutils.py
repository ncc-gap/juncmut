#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""

def juncmut_juncutils(input_file, output_file, cont_list, genome_id, rbamchr, read_num_thres):
    import subprocess
    import shutil
    import pandas as pd
    import os
    import glob
    from junc_utils.utils import proc_star_junction
    #autochom

    infile = input_file
    tmpfile_list = []
    tmpfile1 = output_file + ".tmp1"
    tmpfile_list.append(tmpfile1)

    ##########
    # is this necessary? (YS)    
    indf1 = pd.read_csv(infile, sep='\t', header=None, dtype=str)
    
    if rbamchr == 'chr':
        indf2 = indf1[(indf1.iloc[:,0]=='chr1')|(indf1.iloc[:,0]=='chr2')|(indf1.iloc[:,0]=='chr3')|(indf1.iloc[:,0]=='chr4')|(indf1.iloc[:,0]=='chr5')|(indf1.iloc[:,0]=='chr6')|(indf1.iloc[:,0]=='chr7')|(indf1.iloc[:,0]=='chr8')|(indf1.iloc[:,0]=='chr9')|(indf1.iloc[:,0]=='chr10')|(indf1.iloc[:,0]=='chr11')|(indf1.iloc[:,0]=='chr12')|(indf1.iloc[:,0]=='chr13')|(indf1.iloc[:,0]=='chr14')|(indf1.iloc[:,0]=='chr15')|(indf1.iloc[:,0]=='chr16')|(indf1.iloc[:,0]=='chr17')|(indf1.iloc[:,0]=='chr18')|(indf1.iloc[:,0]=='chr19')|(indf1.iloc[:,0]=='chr20')|(indf1.iloc[:,0]=='chr21')|(indf1.iloc[:,0]=='chr22')|(indf1.iloc[:,0]=='chrX')]
        indf2.to_csv(tmpfile1, sep='\t', header=False, index=False)
    else:
        indf2 = indf1[(indf1.iloc[:,0]=='1')|(indf1.iloc[:,0]=='2')|(indf1.iloc[:,0]=='3')|(indf1.iloc[:,0]=='4')|(indf1.iloc[:,0]=='5')|(indf1.iloc[:,0]=='6')|(indf1.iloc[:,0]=='7')|(indf1.iloc[:,0]=='8')|(indf1.iloc[:,0]=='9')|(indf1.iloc[:,0]=='10')|(indf1.iloc[:,0]=='11')|(indf1.iloc[:,0]=='12')|(indf1.iloc[:,0]=='13')|(indf1.iloc[:,0]=='14')|(indf1.iloc[:,0]=='15')|(indf1.iloc[:,0]=='16')|(indf1.iloc[:,0]=='17')|(indf1.iloc[:,0]=='18')|(indf1.iloc[:,0]=='19')|(indf1.iloc[:,0]=='20')|(indf1.iloc[:,0]=='21')|(indf1.iloc[:,0]=='22')|(indf1.iloc[:,0]=='X')]
        indf2.to_csv(tmpfile1, sep='\t', header=False, index=False)
    #########

    tmpfile2 = output_file + ".tmp2"
    
    if not cont_list:
        proc_star_junction(tmpfile1, tmpfile2, None, read_num_thres, 10, False, False)
        tmpfile_list.append(tmpfile2)
    else:
        n = 1
        cur_infile, cur_outfile = tmpfile1, tmpfile2 + "_" + str(n)
        for cont in cont_list:
            proc_star_junction(cur_infile, cur_outfile, cont, read_num_thres, 10, False, False) #<--reads>=1 adapt 2pass
            tmpfile_list.append(cur_outfile)
            n = n + 1
            cur_infile, cur_outfile = cur_outfile, tmpfile2 + "_" + str(n)

        shutil.copy(cur_infile, tmpfile2)
        tmpfile_list.append(tmpfile2)
        
    annotate_commands = ["junc_utils", "annotate", tmpfile2, output_file, "--genome_id", genome_id, '--gene_model=gencode']
    subprocess.call(annotate_commands)

    for tfile in tmpfile_list:
        os.remove(tfile)
   
 
if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input", metavar = "prefix", default = None, type = str,
                            help = "prefix") 
        
    parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                            help = "folder name of input files") 
        
    parser.add_argument('--control_file', nargs='*', type = str, default = [],
                            help = "control data created by merge_control (default: %(default)s), reads filter>=1")
        
    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                              help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")
    
    parser.add_argument("--rbam_chr_prefix", choices = ["chr", "none"], default = "none",
                              help = "chr prefix used in your bam (default: %(default)s)")
    
    args = parser.parse_args()

    pr = args.input
    folder = args.folder
    cont_list = args.control_file
    genome_id = args.genome_id
    rbamchr = args.rbam_chr_prefix
    print(pr) 
    juncmut_juncutils(pr, folder, cont_list, genome_id, rbamchr)
    


