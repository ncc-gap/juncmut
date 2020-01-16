#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""

def juncmut_juncutils(pr, folder, cont_list, genome_id, rbamchr):
    import subprocess
    import shutil
    import pandas as pd
    import os
    import glob
    from junc_utils import utils
    #autochom
    infile = './data/%s/%s.SJ.out.tab' %(folder,pr)
    f = './data/%s/alterativeSJ_fil_annot/tmp_%s.SJ.out.tab' %(folder,pr)
    
    indf1 = pd.read_csv(infile, sep='\t', header=None, dtype=str)
    
    if rbamchr == 'chr':
        indf2 = indf1[(indf1.iloc[:,0]=='chr1')|(indf1.iloc[:,0]=='chr2')|(indf1.iloc[:,0]=='chr3')|(indf1.iloc[:,0]=='chr4')|(indf1.iloc[:,0]=='chr5')|(indf1.iloc[:,0]=='chr6')|(indf1.iloc[:,0]=='chr7')|(indf1.iloc[:,0]=='chr8')|(indf1.iloc[:,0]=='chr9')|(indf1.iloc[:,0]=='chr10')|(indf1.iloc[:,0]=='chr11')|(indf1.iloc[:,0]=='chr12')|(indf1.iloc[:,0]=='chr13')|(indf1.iloc[:,0]=='chr14')|(indf1.iloc[:,0]=='chr15')|(indf1.iloc[:,0]=='chr16')|(indf1.iloc[:,0]=='chr17')|(indf1.iloc[:,0]=='chr18')|(indf1.iloc[:,0]=='chr19')|(indf1.iloc[:,0]=='chr20')|(indf1.iloc[:,0]=='chr21')|(indf1.iloc[:,0]=='chr22')|(indf1.iloc[:,0]=='chrX')]
        indf2.to_csv(f, sep='\t', header=False, index=False)
    else:
        indf2 = indf1[(indf1.iloc[:,0]=='1')|(indf1.iloc[:,0]=='2')|(indf1.iloc[:,0]=='3')|(indf1.iloc[:,0]=='4')|(indf1.iloc[:,0]=='5')|(indf1.iloc[:,0]=='6')|(indf1.iloc[:,0]=='7')|(indf1.iloc[:,0]=='8')|(indf1.iloc[:,0]=='9')|(indf1.iloc[:,0]=='10')|(indf1.iloc[:,0]=='11')|(indf1.iloc[:,0]=='12')|(indf1.iloc[:,0]=='13')|(indf1.iloc[:,0]=='14')|(indf1.iloc[:,0]=='15')|(indf1.iloc[:,0]=='16')|(indf1.iloc[:,0]=='17')|(indf1.iloc[:,0]=='18')|(indf1.iloc[:,0]=='19')|(indf1.iloc[:,0]=='20')|(indf1.iloc[:,0]=='21')|(indf1.iloc[:,0]=='22')|(indf1.iloc[:,0]=='X')]
        indf2.to_csv(f, sep='\t', header=False, index=False)
    
    file1 = './data/%s/alterativeSJ_fil_annot/%s.SJ.fil.txt' %(folder,pr)
    file2 = './data/%s/alterativeSJ_fil_annot/%s.SJ.fil.annot.txt' %(folder,pr)
        
    t1 = './data/%s/alterativeSJ_fil_annot/tmp_out_%s' %(folder, pr)
    t2 = './data/%s/alterativeSJ_fil_annot/tmp_in_%s' %(folder, pr)
    
    if not cont_list:
        shutil.copy(f, file1)
    
    else:
        #f = './junction/%s/%s.SJ.out.tab' %(args.folder,pr)
        #cont_list : list of cotrol files arg.control_file
        n=1
        for cont in cont_list:
            out = t1+str(n)+'.txt'
            #junc_utils filter --pooled_control_file *.bed.gz input(./junction/A427.SJ.out.tab) output
            utils.proc_star_junction(f, out, cont, 1, 10, False, False) #<--reads>=1 adapt 2pass
            f = t2+str(n)+'.txt'
            shutil.copy(out, f)
            n=n+1 
        shutil.copy(out, file1)
        
    annotate_commands = ["junc_utils", "annotate", file1, file2, "--genome_id", genome_id, '--gene_model=gencode']
    subprocess.call(annotate_commands)
    
    file_list = glob.glob("./data/%s/alterativeSJ_fil_annot/tmp*" %(folder))
    for file in file_list:
        os.remove(file)
    
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
    


