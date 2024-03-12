#! /usr/bin/env python

from annot_utils.run import *
import argparse
from .version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "annot_utils")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    subparsers = parser.add_subparsers()

    ##########
    # gene 

    gene = subparsers.add_parser("gene",
                                  help = "make gene information file")

    gene.add_argument("output_path", metavar = "gene.bed.gz", default = None, type = str,
                      help = "the path to the output")

    gene.add_argument("--gene_model", choices = ["refseq", "gencode"], default = "refseq",
                       help = "gene model (refGene or ensGene) (default: %(default)s)")

    gene.add_argument("--grc", default = False, action = 'store_true',
                      help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    gene.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                      help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    gene.add_argument("--add_ref_id", default = False, action = 'store_true',
                      help = "add refGene ID to the annotation bed file (default: %(default)s)")

    gene.set_defaults(func = gene_main)

    ##########
    # exon 

    exon = subparsers.add_parser("exon",
                                 help = "make exon information file")

    exon.add_argument("output_path", metavar = "exon.bed.gz", default = None, type = str,
                  help = "the path to the output")

    exon.add_argument("--gene_model", choices = ["refseq", "gencode"], default = "refseq",
                   help = "gene model (refGene or ensGene) (default: %(default)s)")

    exon.add_argument("--grc", default = False, action = 'store_true',
                  help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    exon.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                  help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    exon.add_argument("--add_ref_id", default = False, action = 'store_true',
                  help = "add refGene ID to the annotation bed file (default: %(default)s)")

    exon.set_defaults(func = exon_main)

    ##########
    # coding 

    coding = subparsers.add_parser("coding",
                                   help = "make coding information file")

    coding.add_argument("output_path", metavar = "coding.bed.gz", default = None, type = str,
                        help = "the path to the output")

    coding.add_argument("--gene_model", choices = ["refseq", "gencode"], default = "refseq",
                        help = "gene model (refGene or ensGene) (default: %(default)s)")

    coding.add_argument("--grc", default = False, action = 'store_true',
                        help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    coding.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                        help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    coding.add_argument("--add_ref_id", default = False, action = 'store_true',
                        help = "add refGene ID to the annotation bed file (default: %(default)s)")

    coding.set_defaults(func = coding_main)


    ##########
    # junction

    junction = subparsers.add_parser("junction",
                                     help = "make exon-exon junction information file")

    junction.add_argument("output_path", metavar = "junction.bed.gz", default = None, type = str,
                          help = "the path to the output")

    junction.add_argument("--gene_model", choices = ["refseq", "gencode"], default = "refseq",
                          help = "gene model (refGene or ensGene) (default: %(default)s)")

    junction.add_argument("--grc", default = False, action = 'store_true',
                          help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    junction.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                          help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    junction.add_argument("--add_ref_id", default = False, action = 'store_true',
                          help = "add refGene ID to the annotation bed file (default: %(default)s)")

    junction.set_defaults(func = junction_main)


    ##########
    # boundary

    boundary = subparsers.add_parser("boundary",
                                     help = "make boundary information file")

    boundary.add_argument("output_path", metavar = "boudary.bed.gz", default = None, type = str,
                          help = "the path to the output")

    boundary.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                          help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    boundary.add_argument("--grc", default = False, action = 'store_true',
                          help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    boundary.add_argument("--donor_size", metavar = "donor_size", default = "2,6", type = str,
                          help = "splicing donor site size (exonic region size, intronic region size) (default: %(default)s)")

    boundary.add_argument("--acceptor_size", metavar = "acceptor_size", default = "8,1", type = str,
                          help = "splicing donor site size (intronic region size, exonic region size) (default: %(default)s)")

    boundary.set_defaults(func = boundary_main)


    """
    ##########
    # simple_repeat

    simple_repeat = subparsers.add_parser("simple_repeat", help = "make simple repeat information file")
    
    simple_repeat.add_argument("output_path", metavar = "simple_repeat.bed.gz", default = None, type = str,
                               help = "the path to the output")

    simple_repeat.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                               help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    simple_repeat.add_argument("--grc", default = False, action = 'store_true',
                               help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    simple_repeat.set_defaults(func = simple_repeat_main)

    ##########
    """

    return parser


