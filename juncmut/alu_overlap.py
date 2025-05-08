#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import csv

def alu_overlap(input_file, output_file, alu_bed):

    # collect bed
    mut_dict = {}
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            if csvobj["SS_shift1"] == "NA":
                continue

            mut_key = csvobj["Mut_key"]
            mut_key_chr, mut_key_pos, _, _ = mut_key.split(',')
            mut_key_pos = int(mut_key_pos)
            mut_dict[mut_key] = "%s\t%d\t%d\t%s\t%s\n" % (mut_key_chr, mut_key_pos - 1, mut_key_pos, mut_key, csvobj["SJ_strand"])

    with open(output_file+".tmp1.bed", 'w') as hout_tmp:
        for val in mut_dict.values():
            hout_tmp.write(val)

    with open(output_file+".overlap.bed", 'w') as hout_overlap:
        subprocess.run(["bedtools", "intersect", "-a", output_file+".tmp1.bed", "-b", alu_bed, "-wa", "-wb"], stdout = hout_overlap)
    
    #with open(output_file+".nonoverlap.bed", 'w') as hout_nonoverlap:
    #    subprocess.run(["bedtools", "intersect", "-a", output_file+".tmp1.bed", "-b", alu_bed, "-v"], stdout = hout_nonoverlap)

    overlap = {}
    with open(output_file+".overlap.bed") as hin_overlap:
        for line in hin_overlap:
            F = line.rstrip('\n').split('\t')
            mut_key = F[3]
            gene_dir = F[4]
            alu_dir = F[8]

            alu_direction = "-"
            if gene_dir == "+" and alu_dir == "+":
                alu_direction = "sense"
            elif gene_dir == "-" and alu_dir == "-":
                alu_direction = "sense" 
            elif gene_dir == "+" and alu_dir == "-":
                alu_direction = "antisense"
            elif gene_dir == "-" and alu_dir == "+":
                alu_direction = "antisense"

            overlap[mut_key] = {
                "Alu_Chr": F[5],
                "Alu_start_pos": F[6],
                "Alu_end_pos": F[7],
                "Alu_strand": F[8],
                "repName": F[9],
                "repClass": F[10],
                "repFamily": F[11],
                "Alu_direction": alu_direction
            }

    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter = '\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Alu_Chr", "Alu_start_pos", "Alu_end_pos", "Alu_strand", "repName", "repClass", "repFamily", "Alu_direction"
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:
            mut_key = csvobj["Mut_key"]
            if mut_key in overlap:
                csvobj["Alu_Chr"] = overlap[mut_key]["Alu_Chr"]
                csvobj["Alu_start_pos"] = overlap[mut_key]["Alu_start_pos"]
                csvobj["Alu_end_pos"] = overlap[mut_key]["Alu_end_pos"]
                csvobj["Alu_strand"] = overlap[mut_key]["Alu_strand"]
                csvobj["repName"] = overlap[mut_key]["repName"]
                csvobj["repClass"] = overlap[mut_key]["repClass"]
                csvobj["repFamily"] = overlap[mut_key]["repFamily"]
                csvobj["Alu_direction"] = overlap[mut_key]["Alu_direction"]
            else:
                csvobj["Alu_Chr"] = "NA"
                csvobj["Alu_start_pos"] = "NA"
                csvobj["Alu_end_pos"] = "NA"
                csvobj["Alu_strand"] = "NA"
                csvobj["repName"] = "NA"
                csvobj["repClass"] = "NA"
                csvobj["repFamily"] = "NA"
                csvobj["Alu_direction"] = "NA"

            csvwriter.writerow(csvobj)

    os.remove(output_file+".overlap.bed")
    #os.remove(output_file+".nonoverlap.bed")
    os.remove(output_file+".tmp1.bed")

if __name__== "__main__":
    import sys
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    alu_bed = sys.argv[3]
    
    alu_overlap(input_file, output_file, alu_bed)

