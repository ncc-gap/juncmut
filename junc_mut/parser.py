#! /usr/bin/env python

from .run import run_juncm
import argparse
#　

def create_parser():
    
    parser = argparse.ArgumentParser(prog = "juncm") #make a parser

    parser.add_argument("input_file", metavar = "./junction/sample.SJ.out.tab", default = None, type = str,
                        help = "Path to input file") 
    
    parser.add_argument('--control_file', nargs='*', type = str,
                        help = "Path to control data created by merge_control (default: %(default)s)")
    
    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                          help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")
    
    parser.add_argument("--read_num_thres", type = int, default = 3,
                        help = "Remove splicing junctions whose supporting numbers are below this value (default: %(default)s)")
    
    parser.add_argument("--freq_thres", type = int, default = 0.05,
                        help = "Remove splicing junctions whose frequency is below this value (default: %(default)s)")
    
    parser.set_defaults(func = run_juncm) #which def do you use?
    
    
    return parser


