#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import subprocess
import gzip

def define_longest_transcript(splice_type, mut_key_chr, normal_pos, ref_file):
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
            if splice_type == "5'SS+" or splice_type == "3'SS-":
                exonStarts = R[9].split(',')
                if normal_pos in exonStarts:
                    gene2tx_info[ref_tx_id] = mut_key_chr, ref_tx_start, ref_tx_end, ref_gene
            elif splice_type == "5'SS-" or splice_type == "3'SS+":
                exonEnds = R[10].split(',')
                if normal_pos in exonEnds:
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
    genes_list = []
    genes = '---'
    for gene in gene2chr:
        region_list.append("%s:%s-%s" % (gene2chr[gene], gene2start[gene], gene2end[gene]))
        genes_list.append(gene)
        genes = ','.join(genes_list)
        
    return(genes, region_list)

def juncmut_filt_bam_main(input_file, output_file, input_bam, output_bam, genecode_gene_file):

    ex_region_list =[]
    # open a file and make a trnscript list for RNA_Mut True.
    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + ["Gene"])
        csvwriter.writeheader()
        for csvobj in csvreader:
            mut_key_chr = csvobj["Mut_key"].split(',')[0]
            sj_pos = csvobj["SJ_key"].split(':')[1].split('-')
            sj_start = int(sj_pos[0])
            sj_end = int(sj_pos[1])
            
            #o-->
            if "5'SS" in csvobj["Splicing_class"] and csvobj["SJ_strand"] == "+": 
                splice_type = "5'SS+"
                normal_pos = sj_end
            #<--o
            elif "5'SS" in csvobj["Splicing_class"] and csvobj["SJ_strand"] == "-":
                splice_type = "5'SS-"
                normal_pos = sj_start - 1
            #-->o
            elif "3'SS" in csvobj["Splicing_class"] and csvobj["SJ_strand"] == "+":
                splice_type = "3'SS+"
                normal_pos = sj_start - 1
            #o<--
            elif "3'SS" in csvobj["Splicing_class"] and csvobj["SJ_strand"] == "-":
                splice_type = "3'SS-"
                normal_pos = sj_end

            genes, region_list = define_longest_transcript(splice_type, mut_key_chr, str(normal_pos), genecode_gene_file)

            csvobj["Gene"] = genes
            csvwriter.writerow(csvobj)

            ex_region_list.extend(region_list)

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
    output_file = sys.argv[2]
    input_bam = sys.argv[3]
    output_bam = sys.argv[4]
    gencode_gene_file = sys.argv[5]

    juncmut_filt_bam_main(input_file, output_file, input_bam, output_bam, gencode_gene_file)
    
