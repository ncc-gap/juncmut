#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:30:10 2019

@author: genome
"""
def juncmut_annotgnomadsnp(pr, folder, genome_id):
    import pysam
      
    file = '.data/%s/alterativeSJ_mutprediction/%s.SJ.fil.annot.assadjunifreqT.pmut.SJinSJ.txt' %(folder,pr) 
    
    out_file = '.data/%s/alterativeSJ_mutprediction/%s.SJ.fil.annot.assadjunifreqT.pmut.SJinSJ.snp.txt' %(folder,pr)
    
    
    if genome_id == "hg19":
        
        db = "./reference/gnomad.genomes.r2.1.1.sites.vcf.bgz" #chr_prefix is "none".
        tb = pysam.TabixFile(db)
       
        with open(file, 'r') as hin:
            with open(out_file,'w') as hout:
                for line in hin: # one SJ
                    line = line.rstrip('\n')
                    junc_record = line
                    F = line.split('\t')
                    c = (F[0].replace('chr', '')) 
        
                    if c in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]:
                        if F[13] == "-":
                            out_record = line + "\t-\t-\t-\t-\t0\n" 
                            hout.write(out_record)
                        else:
                            P = F[13].split(',') # multiple position per SJ.
                            for i in P:
                                pb = i.split(':')
    
                                for a in pb[2]:
                                    #chr = "chr" + str(c)
                                    chr = str(c) #
                                    rows =tb.fetch(chr, int(pb[0])-1, int(pb[0]))
                                    if not list(rows):
                                        out_record = junc_record + "\t" + str(pb[0]) + "\t" + str(pb[1]) + "\t" + str(a) + "\t" + "-" + "\t0.0\n" 
                                        hout.write(out_record)
        
                                for a in pb[2]:
                                    out_record = ""
                                    #chr = "chr" + str(c)
                                    chr = str(c) #
                                    rows =tb.fetch(chr, int(pb[0])-1, int(pb[0]))
                                    for row in rows:
                                        srow=str(row)
                                        record = srow.split('\t')
                                        
                                        if pb[1] == record[3] and a == record[4]:
                                            allele = record[3]+">"+record[4]
                                            infos = record[7].split(';')
                                            for info in infos:
                                                if info.startswith("AF="):
                                                    freq = float(info.replace("AF=", ''))
        
                                            out_record = junc_record + "\t" + str(pb[0]) + "\t" + str(pb[1]) + "\t" + str(a) + "\t" + allele + "\t" + str(freq) +"\n" 
                                            print(freq)
                                            break                                        
               
                                        out_record = junc_record + "\t" + str(pb[0]) + "\t" + str(pb[1]) + "\t" + str(a) + "\t" + "-" + "\t0.0\n"
                                    hout.write(out_record)
        
                                        
                    elif c in ["Y"]:
                        out_record = junc_record + "\t-\t-\t-\tna\t0\n" 
                        hout.write(out_record)
                    else:
                        out_record = junc_record + "\tPOS\tREF\tMUT_prediction\tsnp\tsnp_freq\n" 
                        hout.write(out_record)
                        
    elif genome_id == "hg38":
        
        db = "./reference/gnomad.genomes.r3.0.sites.vcf.bgz"  #chr_prefix is "chr".
        tb = pysam.TabixFile(db)
       
        with open(file, 'r') as hin:
            with open(out_file,'w') as hout:
                for line in hin: # one SJ
                    line = line.rstrip('\n')
                    junc_record = line
                    F = line.split('\t')
                    c = (F[0].replace('chr', '')) 
        
                    if c in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]:
                        if F[13] == "-":
                            out_record = line + "\t-\t-\t-\t-\t0\n" 
                            hout.write(out_record)
                        else:
                            P = F[13].split(',') # multiple position per SJ.
                            for i in P:
                                pb = i.split(':')
    
                                for a in pb[2]:
                                    chr = "chr" + str(c)
                                    #chr = str(c) #
                                    rows =tb.fetch(chr, int(pb[0])-1, int(pb[0]))
                                    if not list(rows):
                                        out_record = junc_record + "\t" + str(pb[0]) + "\t" + str(pb[1]) + "\t" + str(a) + "\t" + "-" + "\t0.0\n" 
                                        hout.write(out_record)
        
                                for a in pb[2]:
                                    out_record = ""
                                    chr = "chr" + str(c)
                                    #chr = str(c) #
                                    rows =tb.fetch(chr, int(pb[0])-1, int(pb[0]))
                                    for row in rows:
                                        srow=str(row)
                                        record = srow.split('\t')
                                        
                                        if pb[1] == record[3] and a == record[4]:
                                            allele = record[3]+">"+record[4]
                                            infos = record[7].split(';')
                                            for info in infos:
                                                if info.startswith("AF="):
                                                    freq = float(info.replace("AF=", ''))
        
                                            out_record = junc_record + "\t" + str(pb[0]) + "\t" + str(pb[1]) + "\t" + str(a) + "\t" + allele + "\t" + str(freq) +"\n" 
                                            print(freq)
                                            break                                        
               
                                        out_record = junc_record + "\t" + str(pb[0]) + "\t" + str(pb[1]) + "\t" + str(a) + "\t" + "-" + "\t0.0\n"
                                    hout.write(out_record)
        
                                        
                    elif c in ["Y"]:
                        out_record = junc_record + "\t-\t-\t-\tna\t0\n" 
                        hout.write(out_record)
                    else:
                        out_record = junc_record + "\tPOS\tREF\tMUT_prediction\tsnp\tsnp_freq\n" 
                        hout.write(out_record)

if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input", metavar = "prefix", default = None, type = str,
                            help = "Prefix of input file") 
        
    parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                            help = "folder name") 
    
    parser.add_argument("--genome_id", choices = ["hg19", "hg38"], default = "hg19",
                              help = "Genome id used for selecting snp data (default: %(default)s)") 
        
    args = parser.parse_args()
    
    pr = args.input
    folder = args.folder
    genome_id = args.genome_id
    
    juncmut_annotgnomadsnp(pr, folder, genome_id)
