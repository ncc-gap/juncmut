#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida

proc_star_junction(input_file, output_file, control_file, read_num_thres, overhang_thres, remove_annotated, convert_map_splice2):
convert_map_splice2: start position -1 , end position +1
remove_annotated == True to remove a line with column[5] != 0.
"""

def juncmut_juncutils(input_file, output_file, cont_list, genome_id, read_num_thres):
    import subprocess
    import shutil
    import os
    from junc_utils.utils import proc_star_junction
    #autochom


    tmpfile_list = []
    tmpfile1 = output_file + ".tmp1"
    tmpfile_list.append(tmpfile1)

    target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

    with open(input_file, 'r') as hin, open(tmpfile1, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] in target_rnames + ["chr" + x for x in target_rnames]:

                print('\t'.join(F), file = hout)


    tmpfile2 = output_file + ".tmp2"
    
    if not cont_list:
        proc_star_junction(tmpfile1, tmpfile2, None, read_num_thres, 10, False, False)
        tmpfile_list.append(tmpfile2)
    else:
        n = 1
        cur_infile, cur_outfile = tmpfile1, tmpfile2 + "_" + str(n)
        for cont in cont_list:
            proc_star_junction(cur_infile, cur_outfile, cont, read_num_thres, 10, False, False)
            tmpfile_list.append(cur_outfile)
            n = n + 1
            cur_infile, cur_outfile = cur_outfile, tmpfile2 + "_" + str(n)

        shutil.copy(cur_infile, tmpfile2)
        tmpfile_list.append(tmpfile2)
    
    output_path = output_file
    annotate_commands = ["junc_utils", "annotate", tmpfile2, output_path, "--genome_id", genome_id, '--gene_model=gencode']
    subprocess.call(annotate_commands)

    for tfile in tmpfile_list:
        os.remove(tfile)
   
 
if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--input_file", metavar = "input_file", default = None, type = str,
                            help = "prefix") 
        
    parser.add_argument("--output_file", metavar = "output_file", default = "my_sample", type = str,
                            help = "folder name of input files") 
        
    parser.add_argument('--control_file', nargs='*', type = str, default = [],
                            help = "control data created by merge_control (default: %(default)s), reads filter>=1")
        
    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg38",
                              help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")
    
    parser.add_argument("--read_num_thres", metavar = "read_num_thres", default = 1, type = int,
                              help = "chr prefix used in your bam (default: %(default)s)")
    
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    cont_list = args.control_file
    genome_id = args.genome_id
    read_num_thres = args.read_num_thres
    print(input_file)
    
    juncmut_juncutils(input_file, output_file, cont_list, genome_id, read_num_thres)


