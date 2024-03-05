#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_intersect(input_file, output_file, original_sj_file):
    import subprocess
    import os
    import csv

    tmpfile1 = output_file + ".tmp1"
    tmpfile2 = output_file + ".tmp2"

    header = []
    with open(input_file) as hin, open(tmpfile1, 'w') as hout:
        csvreader = csv.reader(hin, delimiter='\t')
        header = next(csvreader)
        for row in csvreader:
            (sj_key_chr, sj_key_start, sj_key_end) = row[header.index("SJ_key")].split(",")
            print("%s\t%d\t%d\t%s" % (sj_key_chr, int(sj_key_start) - 5, int(sj_key_end) + 5, "\t".join(row)), file = hout)

    with open(tmpfile2, 'w') as hout:
        subprocess.run(["bedtools", "intersect", "-a", tmpfile1, "-b", original_sj_file, "-c"], stdout = hout)

    with open(tmpfile2) as hin, open(output_file, 'w') as hout:
        print('\t'.join(header + ["SJ_overlap_count"]), file = hout)
        csvreader = csv.reader(hin, delimiter='\t')
        for row in csvreader:
            print('\t'.join(row[3:]), file = hout)

    os.remove(tmpfile2)
    os.remove(tmpfile1)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    original_sj_file = sys.argv[3]
    juncmut_intersect(input_file, output_file, original_sj_file)
