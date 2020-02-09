#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:39:58 2019

@author: genome
"""
from .juncmut_env import juncmut_env
from .juncmut_juncutils import juncmut_juncutils 
from .juncmut_assadj import juncmut_assadj
from .juncmut_freq import juncmut_freq
from .juncmut_mutpre import juncmut_mutpre
from .juncmut_intersect import juncmut_intersect
from .juncmut_annotgnomadsnp import juncmut_annotgnomadsnp
from .juncmut_annotrnamut import juncmut_annotrnamut

def run_juncmut(args):
    
    # from . import juncmut_env, juncmut_juncutils, juncmut_assadj, juncmut_freq, juncmut_mutpre, juncmut_intersect, juncmut_annotgnomadsnp, juncmut_annotrnamut
    
    # pr = args.input
    
    # otuput_dir = args.output_dir
    
    cont_list = args.control_file
    
    genome_id = args.genome_id
    
    read_num_thres = args.read_num_thres
    
    freq_thres = args.freq_thres
    
    rbamchr = args.rbam_chr_prefix
    
    rbam = args.rbam

     
    juncmut_env(args.output_dir)
   
    juncmut_juncutils(args.input_SJ, args.output_dir, cont_list, genome_id, rbamchr, args.read_num_thres)
   
    juncmut_assadj(args.input_SJ, args.output_dir)
    
    juncmut_freq(args.input_SJ, args.output_dir, read_num_thres, freq_thres)

    """
    juncmut_mutpre(pr, folder, genome_id)

    juncmut_intersect(pr, folder)

    juncmut_annotgnomadsnp(pr, folder, genome_id)
    
    juncmut_annotrnamut(pr, folder, genome_id, rbamchr, rbam)
    """

