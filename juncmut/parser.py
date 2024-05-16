#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from . import __version__
from .run_juncmut import get_main 
from .run_juncmut import annot_main 
from .run_juncmut import filt_bam_main 
from .run_juncmut import filt_main 
from .run_juncmut import sjclass_main 
#from .run_juncmut import validate_main 
import argparse

def create_parser():
    parser = argparse.ArgumentParser(prog = "juncmut")
    parser.add_argument("--version", action = "version", version = __version__)

    subparsers = parser.add_subparsers()
    
    #####
    # get
    get = subparsers.add_parser("get", help = "Identify the mutation to create the abnormal alternative splicing.")
    
    get.add_argument("sj_file", metavar = "sj_file", type = str, help = "Input file(such as sample.SJ.out.tab file generated by STAR)")
    get.add_argument("bam_file", metavar = "bam_file", type = str, help = "Path to RNA bam file")
    get.add_argument("output_file", metavar = "output_file", type = str, help = "Output file")
    get.add_argument("reference", metavar = "reference", type = str, help = "Path to reference genome")
    get.add_argument('--control_file', nargs='*', type = str, help = "Path to control data created by merge_control (default: %(default)s)")
    get.add_argument("--genecode_gene_file", metavar = "genecode_gene_file", required = True, type = str, help = "Path to genecode gene file")
    get.add_argument("--output_bam", metavar = "output_bam", required = True, type = str, help = "Output bam")
    get.add_argument("--read_num_thres", type = int, default = 3, help = "Splicing junctions with reads >= read_num_thres is saved (default: %(default)s)")
    get.add_argument("--freq_thres", type = float, default = 0.05, help = "Splicing junctions with reads >= freq_thres is saved (default: %(default)s)")
    get.add_argument("--mut_num_thres", type = int, default = 1, help = "A mutation with mutation alleles >= mut_num_thres is a true candidate (default: %(default)s)")
    get.add_argument("--mut_freq_thres", type = float, default = 0.05, help = "A mutation with frequency >= mut_freq_thres is a true candidate (default: %(default)s)")
    get.add_argument("--debug", action='store_true', help = "True keeps the intermediate files.")
    
    get.set_defaults(func = get_main)

    #####
    # filt_bam
    filt_bam = subparsers.add_parser("filt_bam", help = "Filt bam by juncmut.txt.")
    
    filt_bam.add_argument("input_file", metavar = "input_file", type = str, help = "juncmut.txt")
    filt_bam.add_argument("bam_file", metavar = "bam_file", type = str, help = "Path to RNA bam file")
    filt_bam.add_argument("output_file", metavar = "output_file", type = str, help = "Output bam file")

    filt_bam.set_defaults(func = filt_bam_main)

    #####
    # filt
    filt = subparsers.add_parser("filt", help = "Filt mutation by SJ_Overlap_Count and gnomAD_AF.")
    
    filt.add_argument("input_file", metavar = "input_file", required = True, type = str, help = "Input file(such as juncmut.txt generated by juncmut)")
    filt.add_argument("output_file", metavar = "output_file", required = True, type = str, help = "Output file")

    filt.set_defaults(func = filt_main)

    #####
    # annot
    annot = subparsers.add_parser("annot", help = "Annotate gnomAD information onto the retrieved mutation.")
    
    annot.add_argument("input_file", metavar = "input_file", type = str, help = "Input file(such as juncmut.txt generated by juncmut)")
    annot.add_argument("output_file", metavar = "output_file", type = str, help = "Output file")
    annot.add_argument("reference", metavar = "reference", type = str, help = "Path to reference genome")
    annot.add_argument("--gnomad", metavar = "gnomad", default = "", type = str, help = "Path to gnomad vcf file")
    annot.add_argument("--debug", action='store_true', help = "True keeps the intermediate files.")
    annot.add_argument("--cgc_file", metavar = "cgc_file", default = "", type = str, help = "Path to cancer gene census file")
    annot.add_argument("--clinvar_file", metavar = "clinvar_file", default = "", type = str, help = "Path to clinvar file")
    annot.add_argument("--acmg_file", metavar = "acmg_file", default = "", type = str, help = "Path to acmg file")
    annot.add_argument("--clinvar_star234_file", metavar = "clinvar_star234_file", default = "", type = str, help = "Path to clinvar_star234 file")
    annot.add_argument("--pancan_file", metavar = "pancan_file", default = "", type = str, help = "Path to pancan file")
    annot.add_argument("--dosage_sensitivity_file", metavar = "dosage_sensitivity_file", default = "", type = str, help = "Path to dosage sensitivity file")
    annot.add_argument("--cgd_file", metavar = "cgd_file", default = "", type = str, help = "Path to cgd file")
    
    annot.set_defaults(func = annot_main)

    #####
    # sjclass
    sjclass = subparsers.add_parser("sjclass", help = "Identify the mutation to create the abnormal alternative splicing.")
    sjclass.add_argument("input_file", metavar = "input_file", type = str, help = "Input file(such as file generated by juncmut")
    sjclass.add_argument("output_file", metavar = "output_file", type = str, help = "Output file")
    sjclass.add_argument("bam_file", metavar = "bam_file", type = str, help = "path to bam.")
    sjclass.add_argument("sj_file", metavar = "sj_file", type = str, help = "path to SJ.out.tab.gz")
    sjclass.add_argument("reference", metavar = "reference", type = str, help = "Path to reference genome")
    sjclass.add_argument("--gencode", metavar = "gencode", default = None, type = str, help = "gencode file")
    sjclass.add_argument("--mane", metavar = "mane", default = None, type = str, help = "mane.gff.transcript_tag.json file")
    sjclass.add_argument("--depth_th", metavar = "depth_th", default = 1, type = int, help = "depth in the intron for classification.")
    sjclass.add_argument("--debug", help = "keep temporary files.", action='store_true')

    sjclass.set_defaults(func = sjclass_main)

    #####
    # validate
    #validate = subparsers.add_parser("validate", help = "Validate the mutations on DNA.")
    #
    #validate.add_argument("-input_file", metavar = "input_file", required = True, type = str, help = "Input file obtained by the juncmut get")
    #validate.add_argument("-output_file", metavar = "output_file", required = True, type = str, help = "Output file")
    #validate.add_argument("-reference", metavar = "reference", required = True, type = str, help = "Path to reference genome")
    #validate.add_argument("-dna_bam", metavar = "dna_bam", required = True, type = str, help = "Path to DNA bam file")
    #validate.add_argument("-mut_num_thres", type = int, default = 1, help = "A mutation with mutation alleles >= mut_num_thres is a true candidate (default: %(default)s)")
    #validate.add_argument("-mut_freq_thres", type = float, default = 0.05, help = "A mutation with frequency >= mut_freq_thres is a true candidate (default: %(default)s)")
    #
    #validate.set_defaults(func = validate_main)
    
    #####
    
    return parser
