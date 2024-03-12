#! /usr/bin/env python

from __future__ import print_function
import sys, pkg_resources 


def make_ucsc2grc(genome_id):

    if genome_id == "hg19":
        sequence_report = pkg_resources.resource_filename("annot_utils", "data/hg19/GCF_000001405.13.assembly.txt")
    elif genome_id == "hg38":
        sequence_report = pkg_resources.resource_filename("annot_utils", "data/hg38/GCF_000001405.26.assembly.txt")
    elif genome_id == "mm10":
        sequence_report = pkg_resources.resource_filename("annot_utils", "data/mm10/GCF_000001635.20.assembly.txt")
    else:
        print("genome_id shoud be hg19, hg38 or mm10", file = sys.stderr)
        sys.exit(1)


    ucsc2grch = {}
    with open(sequence_report, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            if line == "": continue
            F = line.rstrip('\n\r').split('\t')
            if F[9] == "na": continue

            if F[4].startswith('CM'):
                ucsc2grch[F[9]] = F[2]
            else:
                ucsc2grch[F[9]] = F[4]

    ucsc2grch["chrM"] = "MT"

    return ucsc2grch

