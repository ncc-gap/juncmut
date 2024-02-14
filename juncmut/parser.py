#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .run_juncmut import get_main 
from .run_juncmut import annot_main 
from .run_juncmut import validate_main 
import argparse

def create_parser():
    
    parser = argparse.ArgumentParser(prog = "juncmut")
    subparsers = parser.add_subparsers()
    
    #####
    # get
    get = subparsers.add_parser("get", help = "Identify the mutation to create the abnormal alternative splicing.")
    
    get.add_argument("-input_file", metavar = "input_file", required = True, type = str, help = "Input file(such as sample.SJ.out.tab file generated by STAR)")
    get.add_argument("-output_file", metavar = "output_file", required = True, type = str, help = "Output file")
    get.add_argument("-reference", metavar = "reference", required = True, type = str, help = "Path to reference genome")
    get.add_argument("-rna_bam", metavar = "rna_bam", required = True, type = str, help = "Path to RNA bam file")
    get.add_argument('-control_file', nargs='*', type = str, help = "Path to control data created by merge_control (default: %(default)s)")
    get.add_argument("-read_num_thres", type = int, default = 3, help = "Splicing junctions with reads >= read_num_thres is saved (default: %(default)s)")
    get.add_argument("-freq_thres", type = float, default = 0.05, help = "Splicing junctions with reads >= freq_thres is saved (default: %(default)s)")
    get.add_argument("-mut_num_thres", type = int, default = 1, help = "A mutation with mutation alleles >= mut_num_thres is a true candidate (default: %(default)s)")
    get.add_argument("-mut_freq_thres", type = float, default = 0.05, help = "A mutation with frequency >= mut_freq_thres is a true candidate (default: %(default)s)")
    get.add_argument("-genecode_gene_file", metavar = "genecode_gene_file", required = True, type = str, help = "Path to genecode gene file")
    get.add_argument("-output_bam", metavar = "output_bam", required = True, type = str, help = "Output bam")
    get.add_argument("-debug", metavar = "debug", default = "False", type = str, help = "True keeps the intermediate files.")
    
    get.set_defaults(func = get_main)
    
    #####
    # annot
    annot = subparsers.add_parser("annot", help = "Identify the mutation to create the abnormal alternative splicing.")
    
    annot.add_argument("-input_file", metavar = "input_file", required = True, type = str, help = "Input file(such as juncmut.txt generated by juncmut)")
    annot.add_argument("-output_file", metavar = "output_file", required = True, type = str, help = "Output file")
    annot.add_argument("-reference", metavar = "reference", required = True, type = str, help = "Path to reference genome")
    annot.add_argument("-gnomad", metavar = "gnomad", required = True, type = str, help = "Path to gnomad vcf file")
    annot.add_argument("-debug", metavar = "debug", default = "False", type = str, help = "True keeps the intermediate files.")
    
    annot.set_defaults(func = annot_main)
    
    #####
    # validate
    validate = subparsers.add_parser("validate", help = "Validate the mutations on DNA.")
    
    validate.add_argument("-input_file", metavar = "input_file", required = True, type = str, help = "Input file obtained by the juncmut get")
    validate.add_argument("-output_file", metavar = "output_file", required = True, type = str, help = "Output file")
    validate.add_argument("-reference", metavar = "reference", required = True, type = str, help = "Path to reference genome")
    validate.add_argument("-dna_bam", metavar = "dna_bam", required = True, type = str, help = "Path to DNA bam file")
    validate.add_argument("-mut_num_thres", type = int, default = 1, help = "A mutation with mutation alleles >= mut_num_thres is a true candidate (default: %(default)s)")
    validate.add_argument("-mut_freq_thres", type = float, default = 0.05, help = "A mutation with frequency >= mut_freq_thres is a true candidate (default: %(default)s)")
    
    validate.set_defaults(func = validate_main)
    
    #####
    
    return parser
