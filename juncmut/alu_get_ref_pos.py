#! /usr/bin/env python

import csv
import pysam
from Bio import Align

def get_seq(reference, chr, start, end):

    seq = ""
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        # if item[0] == ">": continue
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')

    return seq

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

def alu_ref_alu_pos(input_file, output_file, ref, alu):
    
    if alu == "AluJ0":
        # 283
        Alu_ref = """ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggaggattgcttgagcc
            caggagttcgagaccagcctgggcaacatagcgagaccccgtctctacaaaaaatacaaaaattagccggg
            cgtggtggcgcgcgcctgtagtcccagctactcgggaggctgaggcaggaggatcgcttgagcccaggagt
            tcgaggctgcagtgagctatgatcgcgccactgcactccagcctgggcgacagagcgagaccctgtctca"""
    elif alu == "AluJ0pA":
        # len = 289
        Alu_ref = """ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggaggattgcttgagcc
            caggagttcgagaccagcctgggcaacatagcgagaccccgtctctacaaaaaatacaaaaattagccggg
            cgtggtggcgcgcgcctgtagtcccagctactcgggaggctgaggcaggaggatcgcttgagcccaggagt
            tcgaggctgcagtgagctatgatcgcgccactgcactccagcctgggcgacagagcgagaccctgtctcaaaaaaa"""
    elif alu == "AluYb8pA":
        Alu_ref = """ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggtggatcatgaggtca
            ggagatcgagaccatcctggctaacaaggtgaaaccccgtctctactaaaaatacaaaaaattagccggg
            cgcggtggcgggcgcctgtagtcccagctactcgggaggctgaggcaggagaatggcgtgaacccgggaa
            gcggagcttgcagtgagccgagattgcgccactgcagtccgcagtccggcctgggcgacagagcgagactccgtctcaaaaaaa"""

    elif alu == "Alu_Funakoshi":
        Alu_ref = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA" + \
            "GGCGGGCGGATCACTTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACAT" + \
            "GGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGG" + \
            "CGCGTGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCT" + \
            "TGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACT" + \
            "CCAGCCTGGGCGACAGAGCGAGACTCTGTCTCAAAAAAAAAAAAAAAAAA"

    #else:
    #    print("No alu parameter.")
   
    Alu_ref = Alu_ref.replace(' ', '').replace('\n', '').upper()

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter = '\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Juncmut_primary_SS_in_Alu", "Juncmut_secondary_SS_in_Alu", "Mut_pos_in_Alu"
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:
            if csvobj["Alu_Chr"] == "NA":
                csvobj["Juncmut_primary_SS_in_Alu"] = "NA"
                csvobj["Juncmut_secondary_SS_in_Alu"] = "NA"
                csvobj["Mut_pos_in_Alu"] = "NA"
                csvwriter.writerow(csvobj)
                continue

            secondary_ss = csvobj["Juncmut_secondary_SS"]
            mut_pos = int(csvobj["Mut_key"].split(',')[1])
            mut_ref = csvobj["Mut_key"].split(',')[2]
            primary_ss = int(csvobj["Juncmut_primary_SS"])
            
            alu_chr = csvobj["Alu_Chr"]
            alu_start_pos = int(csvobj["Alu_start_pos"]) + 1
            alu_end_pos = int(csvobj["Alu_end_pos"])
            alu_strand = csvobj["Alu_strand"]

            csvobj["Juncmut_primary_SS_in_Alu"] = "---"
            csvobj["Juncmut_secondary_SS_in_Alu"] = "---"
            csvobj["Mut_pos_in_Alu"] = "---"

            if csvobj["repFamily"] != "Alu" or primary_ss < alu_start_pos or primary_ss > alu_end_pos:
                csvwriter.writerow(csvobj)
                continue

            if alu_strand == '+':
                rel_primary_ss = primary_ss - alu_start_pos
                rel_mut_pos = mut_pos - alu_start_pos
                rel_secondary_ss = int(secondary_ss) - alu_start_pos if secondary_ss != 'NA' else "---"
            else:
                rel_primary_ss = alu_end_pos - primary_ss
                rel_mut_pos = alu_end_pos - mut_pos
                rel_secondary_ss = alu_end_pos - int(secondary_ss) if secondary_ss != 'NA' else "---"

            rep_seq = get_seq(ref, alu_chr, alu_start_pos, alu_end_pos)
            
            if alu_strand == '-': 
                rep_seq = reverse_complement(rep_seq)
                mut_ref = reverse_complement(mut_ref)

            alignments = aligner.align(rep_seq, Alu_ref)
            
            query_seq_pos = alignments[0].coordinates[0, 0] -1 
            reference_seq_pos = alignments[0].coordinates[1, 0] -1

            for i in range(len(alignments[0][0])):
                if alignments[0][0][i] != "-": query_seq_pos = query_seq_pos + 1
                if alignments[0][1][i] != "-": reference_seq_pos = reference_seq_pos + 1
                if query_seq_pos == rel_mut_pos and alignments[0][0][i] != "-": 
                    csvobj["Mut_pos_in_Alu"] = reference_seq_pos + 1
                    if alignments[0][0][i] != mut_ref:
                        raise Exception(f'Sequence inconsistent, {csvobj["Mut_key"]}, {i}')
                if query_seq_pos == rel_primary_ss and alignments[0][0][i] != "-":
                    csvobj["Juncmut_primary_SS_in_Alu"] = reference_seq_pos + 1
                if query_seq_pos == rel_secondary_ss and alignments[0][0][i] != "-":
                    csvobj["Juncmut_secondary_SS_in_Alu"] = reference_seq_pos + 1

            #if mod_primary_ss == "---" or mod_mut_pos == "---":
            #    continue
            csvwriter.writerow(csvobj)

if __name__== "__main__":
    import sys
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ref = sys.argv[3]
    alu = sys.argv[4]
    
    alu_ref_alu_pos(input_file, output_file, ref, alu)
