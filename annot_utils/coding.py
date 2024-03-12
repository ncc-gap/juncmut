#! /usr/bin/env python

from __future__ import print_function
import sys, gzip, subprocess, pkg_resources 
from . import chr_name, utils

def make_coding_info(output_file, gene_model, genome_id, is_grc, add_ref_id):

    # create UCSC to GRC chr name corresponding table
    ucsc2grc = {} 
    if is_grc:
        ucsc2grc = chr_name.make_ucsc2grc(genome_id)

    ucsc_gene_file = utils.set_ucsc_gene_file(genome_id, gene_model)


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
            exon_starts = F[9].split(',')
            exon_ends = F[10].split(',')
            strand = F[3]
            exonNum = int(F[8])
            gene = F[1]
            symbol = F[12]

            size = 0
            for i in range(len(exon_starts) - 1):
                size = size + int(exon_ends[i]) - int(exon_starts[i])


            gene_print_name = "---"
            if gene_model == "refseq":
                if add_ref_id:
                    gene_print_name = symbol + '(' + gene_id + ')'
                else:
                    gene_print_name = symbol
            elif gene_model == "gencode":
                gene_print_name = gene_id


            # for coding genes
            if F[6] != F[7]: 
                for i in range(len(exon_starts) - 1):
                    if i >= 1:
                        key = chr + '\t' + exon_ends[i - 1] + '\t' + exon_starts[i] + '\t' + gene_print_name
                        print(key + '\t' + "intron" + '\t' + strand, file = hout)
    
                    if min(int(exon_ends[i]), int(F[6])) - int(exon_starts[i]) > 0:
                        key = chr + '\t' + exon_starts[i] + '\t' + str(min(int(exon_ends[i]), int(F[6]))) + '\t' + gene_print_name
                        if strand == '+':
                            print(key + '\t' + "5UTR" + '\t' + '+', file = hout)
                        else:
                            print(key + '\t' + "3UTR" + '\t' + '-', file = hout)

                    if min(int(exon_ends[i]), int(F[7])) - max(int(exon_starts[i]), int(F[6])) > 0:
                        key = chr + '\t' + str(max(int(exon_starts[i]), int(F[6]))) + '\t' + str(min(int(exon_ends[i]), int(F[7]))) + '\t' + gene_print_name
                        print(key + '\t' + "coding" + '\t' + strand, file = hout)

                    if int(exon_ends[i]) - max(int(F[7]), int(exon_starts[i])) > 0:
                        key = chr + '\t' + str(max(int(F[7]), int(exon_starts[i]))) + '\t' + exon_ends[i] + '\t' + gene_print_name
                        if strand == '+':
                            print(key + '\t' + "3UTR" + '\t' + '+', file = hout)
                        else:
                            print(key + '\t' + "5UTR" + '\t' + '-', file = hout)

            # for noncoding genes
            if F[6] == F[7]: 
                for i in range(len(exon_starts) - 1):
                    if i >= 1:
                        key = chr + '\t' + exon_ends[i - 1] + '\t' + exon_starts[i] + '\t' + gene_print_name
                        print(key + '\t' + "intron" + '\t' + strand, file = hout)

                    key = chr + '\t' + exon_starts[i] + '\t' + exon_ends[i] + '\t' + gene_print_name
                    print(key + '\t' + "noncoding" + '\t' + strand, file = hout)

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


