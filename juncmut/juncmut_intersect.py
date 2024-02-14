#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_intersect(input_file, output_file, original_sj_file):
    import subprocess
    import os

    tmpfile1 = output_file + ".tmp1"
    tmpfile2 = output_file + ".tmp2"

    header = []
    with open(input_file) as hin, open(tmpfile1, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if header == []:
                header.extend(F)
                continue
            chrom = F[0]
            start = int(F[1])
            end = int(F[2])
            print("%s\t%d\t%d\t%s" % (chrom, start - 5, end + 5, "\t".join(F)), file = hout)

    with open(tmpfile2, 'w') as hout:
        subprocess.run(["bedtools", "intersect", "-a", tmpfile1, "-b", original_sj_file, "-c"], stdout = hout)

    with open(tmpfile2) as hin, open(output_file, 'w') as hout:
        print('\t'.join(header + ["SJ_overlap_count"]), file = hout)
        for line in hin:
            F = line.rstrip('\n').split('\t')
            print('\t'.join(F[3:]), file = hout)

    os.remove(tmpfile2)
    os.remove(tmpfile1)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    original_sj_file = sys.argv[3]
    juncmut_intersect(input_file, output_file, original_sj_file)
