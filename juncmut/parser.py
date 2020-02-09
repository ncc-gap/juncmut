#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:30:39 2019

@author: genome
"""

from .run_juncmut import run_juncmut
import argparse
#　

def create_parser():
    
    parser = argparse.ArgumentParser(prog = "juncmut") #make a parser

    parser.add_argument("input", metavar = "./sampleFolder/sample.SJ.out.tab", default = None, type = str,
                        help = "Prefix of input file") 
    
    parser.add_argument("output_dir", metavar = "output_dir", default = None, type = str,
                        help = "Path to output directory") 
    
    parser.add_argument('--control_file', nargs='*', type = str,
                        help = "Path to control data created by merge_control (default: %(default)s)")
    
    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                          help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")
    
    parser.add_argument("--read_num_thres", type = int, default = 3,
                        help = "Remove splicing junctions whose supporting numbers are below this value (default: %(default)s)")
    
    parser.add_argument("--freq_thres", type = float, default = 0.05,
                        help = "Remove splicing junctions whose frequency is below this value (default: %(default)s)")
    
    parser.add_argument("--rbam_chr_prefix", choices = ["chr", "none"], default = "none",
                              help = "chr prefix used in your bam (default: %(default)s)")
    
    parser.add_argument("--rbam", metavar = "RNAseq_bam_list", default = None, type = str,
                            help = "A file:list of Path to RNAseq bam folder.") 
    
    parser.set_defaults(func = run_juncmut) #which def do you use?
    
    
    
    return parser
