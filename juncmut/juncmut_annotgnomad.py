#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def juncmut_annotgnomad(input_file, output_file, gnomad, genome_id):
    import csv
    import os
    import pysam

    db = gnomad
    tb = pysam.TabixFile(db)

    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "gnomAD", "gnomAD_AF"
        ])
        csvwriter.writeheader()
        for csvobj in csvreader: # one SJ
            (mut_chr, mut_pos, mut_ref, mut_alt) = csvobj["Mut_key"].split(',')
            mut_pos = int(mut_pos)

            remove_chr = mut_chr.replace('chr', '')
            if remove_chr == "Y":
                csvobj["gnomAD"] = "na"
                csvobj["gnomAD_AF"] = 0
                csvwriter.writerow(csvobj)
                continue

            if genome_id == "hg19":
                chr = remove_chr
            else:
                chr = "chr" + remove_chr

            #skip val
            rows = tb.fetch(chr, mut_pos - 1, mut_pos)

            cur_AF = 0.0
            cur_allele = "-"
            if rows is None:
                rows = []

            for row in rows:
                record = row.split('\t')
                if mut_ref == record[3] and mut_alt == record[4]:
                    allele = "%s>%s" % (record[3], record[4])
                    for info in record[7].split(';'):
                        if info.startswith("AF="):
                            cur_AF = float(info.replace("AF=", ''))
                            cur_allele = allele
                    break

            csvobj["gnomAD"] = cur_allele
            csvobj["gnomAD_AF"] = cur_AF
            csvwriter.writerow(csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gnomad = sys.argv[3]

    juncmut_annotgnomad(input_file, output_file, gnomad, "hg38")
