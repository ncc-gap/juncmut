#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import subprocess
import gzip

def define_longest_transcript(splice_type, mut_key_chr, normal_pos):
    gene2tx_info = {}
    with gzip.open(ref_file, 'rt') as hin:
        for record in hin:
            R = record.rstrip('\n').split('\t')
            ref_chr = R[2]
            if mut_key_chr != ref_chr:
                continue
            ref_tx_id = R[1]
            ref_tx_start = R[4]
            ref_tx_end = R[5]
            ref_gene = R[12]
            if splice_type == "Donor+" or splice_type == "Acceptor-":
                exonStarts = R[9].split(',')
                if str(normal_pos) in exonStarts:
                    gene2tx_info[ref_tx_id] = mut_key_chr, ref_tx_start, ref_tx_end, ref_gene
            elif splice_type == "Donor-" or splice_type == "Acceptor+":
                exonEnds = R[10].split(',')
                if str(normal_pos) in exonEnds:
                    gene2tx_info[ref_tx_id] = mut_key_chr, ref_tx_start, ref_tx_end, ref_gene

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

def juncmut_filt_bam(input_file, input_bam, output_bam):

    ex_region_list =[]
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            mutkey_chr = csvobj["Mut_key"].split(',')[0]
            sjkey_pos = csvobj["SJ_key"].split(':')[1].split('-')
            sjkey_start = int(sj_pos[0])
            sjkey_end = int(sj_pos[1])

            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]
            if splice_type == "Donor+" or splice_type == "Acceptor-": 
                normal_pos = sjkey_end
            elif splice_type == "Donor-" or splice_type == "Acceptor+":
                normal_pos = sjkey_start - 1

            ex_region_list.extend(define_longest_transcript(splice_type, mutkey_chr, normal_pos))

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

    juncmut_filt_bam(input_file, input_bam, output_bam)
