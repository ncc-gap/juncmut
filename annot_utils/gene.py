#! /usr/bin/env python

from __future__ import print_function
import sys, gzip, subprocess 
from . import chr_name, utils

def make_gene_info(output_file, gene_model, genome_id, is_grc, add_ref_id, ucsc_gene_file):

    # create UCSC to GRC chr name corresponding table
    ucsc2grc = {}
    if is_grc:
        ucsc2grc = chr_name.make_ucsc2grc(genome_id)

    #ucsc_gene_file = utils.set_ucsc_gene_file(genome_id, gene_model)

    hout = open(output_file + ".unsorted.tmp", 'w')
    with gzip.open(ucsc_gene_file, 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if is_grc:
                if F[2] not in ucsc2grc: continue
                chr = ucsc2grc[F[2]]
            else:
                chr = F[2]
            # chr = ucsc2grc[F[2]] if F[2] in ucsc2grc else F[2]

            gene_id = F[1]
            gene_start = F[4]
            gene_end = F[5]
            strand = F[3]
            symbol = F[12]
            exon_starts = F[9].split(',')
            exon_ends = F[10].split(',')


            gene_print_name = "---"
            if gene_model == "refseq":
                if add_ref_id:
                    gene_print_name = symbol + '(' + gene_id + ')'
                else:
                    gene_print_name = symbol
            elif gene_model == "gencode":
                gene_print_name = gene_id
            else:
                print("the value of gene_model should be refseq or gencode", file = sys.stderr)
                sys.exit(1)

            key = chr + '\t' + gene_start + '\t' + gene_end
            print(key + '\t' + gene_print_name + '\t' + "0" + '\t' + strand, file = hout)

    hout.close()


    hout = open(output_file + ".sorted.tmp", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".unsorted.tmp"], stdout = hout)
    hout.close()

    hout = open(output_file, 'w')
    subprocess.check_call(["bgzip", "-f", "-c", output_file + ".sorted.tmp"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", output_file])


    subprocess.check_call(["rm", "-rf", output_file + ".unsorted.tmp"])
    subprocess.check_call(["rm", "-rf", output_file + ".sorted.tmp"])


