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

    target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

    tmpfile1 = output_file + ".tmp1"
    with open(input_file, 'r') as hin, open(tmpfile1, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] in target_rnames + ["chr" + x for x in target_rnames]:
                print('\t'.join(F), file = hout)

    tmpfile_list = []
    tmpfile_list.append(tmpfile1)

    tmpfile2 = output_file + ".tmp2"
    if not cont_list:
        proc_star_junction(tmpfile1, tmpfile2, None, read_num_thres, 10, False, False)
        tmpfile_list.append(tmpfile2)
    else:
        cur_infile = tmpfile1
        for n,cont in enumerate(cont_list):
            cur_outfile = "%s_%d" % (tmpfile2, n)
            proc_star_junction(cur_infile, cur_outfile, cont, read_num_thres, 10, False, False)
            tmpfile_list.append(cur_outfile)
            cur_infile = cur_outfile

        shutil.copy(cur_infile, tmpfile2)
        tmpfile_list.append(tmpfile2)

    annotate_commands = ["junc_utils", "annotate", tmpfile2, output_file, "--genome_id", genome_id, '--gene_model=gencode']
    subprocess.call(annotate_commands)
    
    for tfile in tmpfile_list:
        os.remove(tfile)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    cont_list = sys.argv[3].split(",")
    
    juncmut_juncutils(input_file, output_file, cont_list, "hg38", 3)
