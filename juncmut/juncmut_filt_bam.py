#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import subprocess
import gzip

def load_gencode(target_chr, target_pos, gencode_gene_file):
    if gencode_gene_file.endswith("gtf.gz"):
        key_value_split = " "
    elif gencode_gene_file.endswith("gff3.gz"):
        key_value_split = "="
    else:
        raise Exception("unexcepted gencode gene file")

    gencode_tx_info = {}
    with gzip.open(gencode_gene_file, 'rt') as hin:
        for line in hin:
            if line.startswith("#"):
                continue
            F = line.rstrip('\n').split('\t')

            chr = F[0]
            feature_type = F[2]
            start = int(F[3])-1
            end = int(F[4])
            strand = F[6]

            if feature_type == "transcript":
                transcript_id = None
                if chr != target_chr or target_pos < start or target_pos > end:
                    continue
            elif feature_type == "exon":
                if transcript_id is None:
                    continue
            else:
                continue

            add_info = {"MANE": "---", }
            for item in F[8].rstrip(";").split(";"):
                (key, value) = item.strip(" ").replace('"', '').split(key_value_split)
                if key in ["transcript_id", "gene_name", "exon_number"]:
                    add_info[key] = value
                elif key == "tag" and "MANE" in value:
                    add_info["MANE"] = value

            if feature_type == "transcript":
                transcript_id = add_info["transcript_id"]
                gencode_tx_info[transcript_id] = {
                    "Chr": chr,
                    "Transcript_start": start,
                    "Transcript_end": end,
                    "Strand": strand,
                    "Gene_name": add_info["gene_name"],
                    "MANE": add_info["MANE"],
                    "Exon_starts": [],
                    "Exon_ends": [],
                }
            else:
                gencode_tx_info[transcript_id]["Exon_starts"].append(start)
                gencode_tx_info[transcript_id]["Exon_ends"].append(end)

    return gencode_tx_info

def define_longest_transcript(splice_type, mut_key_chr, normal_pos, gencode_gene_file):

    gencode_tx_info = load_gencode(mut_key_chr, normal_pos, gencode_gene_file)

    gene2tx_info = {}
    for transcript_id in gencode_tx_info:
        if splice_type == "Donor+" or splice_type == "Acceptor-":
            if not normal_pos in gencode_tx_info[transcript_id]["Exon_starts"]:
                continue
        elif splice_type == "Donor-" or splice_type == "Acceptor+":
            if not normal_pos in gencode_tx_info[transcript_id]["Exon_ends"]:
                continue
        gene2tx_info[transcript_id] = mut_key_chr, gencode_tx_info[transcript_id]["Transcript_start"], gencode_tx_info[transcript_id]["Transcript_end"], gencode_tx_info[transcript_id]["Gene_name"]

    gene2chr = {}
    gene2start = {}
    gene2end = {}
    for id in gene2tx_info:
        tx_chr, tx_start, tx_end, gene = gene2tx_info[id]
        if gene in gene2chr:
            if gene2chr[gene] != tx_chr:
                print('Wrong.')
        else:
            gene2chr[gene] = tx_chr
                    
        if gene in gene2start:
            if int(gene2start[gene]) > int(tx_start):
                gene2start[gene] = tx_start
        else:
            gene2start[gene] = tx_start
        
        if gene in gene2end:
            if int(gene2end[gene]) < int(tx_end):
                gene2end[gene] = tx_end
        else:
            gene2end[gene] = tx_end

    region_list = []
    for gene in gene2chr:
        region_list.append("%s:%s-%s" % (gene2chr[gene], gene2start[gene], gene2end[gene]))

    return region_list

def juncmut_filt_bam(input_file, input_bam, output_bam, genecode_gene_file):

    ex_region_list =[]
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            mutkey_chr = csvobj["Mut_key"].split(',')[0]
            sjkey_pos = csvobj["SJ_key"].split(':')[1].split('-')
            sjkey_start = int(sjkey_pos[0])
            sjkey_end = int(sjkey_pos[1])

            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]
            if splice_type == "Donor+" or splice_type == "Acceptor-": 
                normal_pos = sjkey_end
            elif splice_type == "Donor-" or splice_type == "Acceptor+":
                normal_pos = sjkey_start - 1

            ex_region_list.extend(define_longest_transcript(splice_type, mutkey_chr, normal_pos, genecode_gene_file))

    if ex_region_list: 
        # initialize the file
        open(output_bam + ".tmp.unsorted.sam", 'w').close()
        with open(output_bam + ".tmp.unsorted.sam", 'a') as hout:
            for region in sorted(list(set(ex_region_list))):
                subprocess.check_call(["samtools", "view", input_bam, region], stdout = hout, stderr = subprocess.DEVNULL)
        
        with open(output_bam + ".tmp.unsorted.rmdup.sam", 'w') as hout:
            subprocess.check_call(["sort", "-u", output_bam + ".tmp.unsorted.sam"], stdout = hout)
        
        with open(output_bam + ".tmp.unsorted2.sam", 'w') as hout:
            subprocess.check_call(["samtools", "view", "-H", input_bam], stdout = hout, stderr = subprocess.DEVNULL)
        
        with open(output_bam + ".tmp.unsorted2.sam", 'a') as hout:
            subprocess.check_call(["cat", output_bam + ".tmp.unsorted.rmdup.sam"], stdout = hout)
        
        with open(output_bam + ".tmp.unsorted2.bam", 'w') as hout:
            subprocess.check_call(["samtools", "view", "-hbS", output_bam + ".tmp.unsorted2.sam"], stdout = hout)
        
        with open(output_bam, 'w') as hout:
            subprocess.check_call(["samtools", "sort", output_bam + ".tmp.unsorted2.bam"], stdout = hout)
        
        subprocess.check_call(["samtools", "index", output_bam])
        
        os.remove(output_bam + ".tmp.unsorted.sam")
        os.remove(output_bam + ".tmp.unsorted.rmdup.sam")
        os.remove(output_bam + ".tmp.unsorted2.sam")
        os.remove(output_bam + ".tmp.unsorted2.bam")
    
    else:
        # initialize the file
        with open(output_bam + ".tmp.unsorted.sam", 'w') as hout:
            subprocess.check_call(["samtools", "view", "-H", input_bam], stdout = hout, stderr = subprocess.DEVNULL)
        with open(output_bam, 'w') as hout:
            subprocess.check_call(["samtools", "view", "-hbS", output_bam + ".tmp.unsorted.sam"], stdout = hout)
        subprocess.check_call(["samtools", "index", output_bam])
        os.remove(output_bam + ".tmp.unsorted.sam")

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    input_bam = sys.argv[2]
    output_bam = sys.argv[3]
    genecode_gene_file = sys.argv[4]

    juncmut_filt_bam(input_file, input_bam, output_bam, genecode_gene_file)
