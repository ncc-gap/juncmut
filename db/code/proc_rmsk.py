#! /usr/bin/env python3

import sys, gzip

rmsk_file = sys.argv[1]
output_file = sys.argv[2]

with gzip.open(rmsk_file, 'rb') as hin, open(output_file, "w") as hout:
    next(hin)
    for line in hin:
        F = line.decode().rstrip('\n').split('\t')

        if F[11] == "SINE":
            hout.write('\t'.join([F[5], F[6], F[7], F[9], F[10], F[11], F[12]])+ '\n')

#rmsk_file="rmsk.txt.gz"
#output_file="rmsk_sine.bed"