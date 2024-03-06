#!/usr/bin/env python3

import csv

# the first column name should be "Gene"
def annot_gene_list(input_file, output_file, gene_list_file, column_name):
    gene_list = {}
    with open(gene_list_file, 'r') as hin:
        for F in csv.DictReader(hin, delimiter = '\t'):
            gene_list[F["Gene"]] = 1

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames + [column_name]
        print('\t'.join(header), file = hout)
        for F in dreader:
            is_gene = "TRUE" if F["Gene"] in gene_list else "FALSE"
            print('\t'.join(F.values()) + '\t' + is_gene, file = hout)
 
if __name__== "__main__":
    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gene_list_file = sys.argv[3]
    column_name = sys.argv[4]

    annot_gene_list(input_file, output_file, gene_list_file, column_name)
