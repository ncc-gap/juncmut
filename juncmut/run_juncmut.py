#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:39:58 2019

@author: genome
"""

def run_juncmut(args):
    
    from . import juncmut_env, juncmut_juncutils, juncmut_assadj, juncmut_freq, juncmut_mutpre, juncmut_intersect, juncmut_annotgnomadsnp, juncmut_annotrnamut
    
    pr = args.input
    
    folder = args.folder
    
    cont_list = args.control_file
    
    genome_id = args.genome_id
    
    read_num_thres = args.read_num_thres
    
    freq_thres = args.freq_thres
    
    rbamchr = args.rbam_chr_prefix
    
    rbam = args.rbam
     
    
    juncmut_env(folder)
    
    juncmut_juncutils(pr, folder, cont_list, genome_id)
    
    juncmut_assadj(pr, folder)
    
    juncmut_freq(pr, folder ,read_num_thres, freq_thres)

    juncmut_mutpre(pr, folder, genome_id)

    juncmut_intersect(pr, folder)

    juncmut_annotgnomadsnp(pr, folder, genome_id)
    
    juncmut_annotrnamut(pr, folder, genome_id, rbamchr, rbam)
