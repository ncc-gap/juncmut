#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:30:10 2019

@author: genome
"""
def juncmut_annotgnomadsnp(input_file, output_file, gnomad_path, genome_id):
    import pysam
      
    db = gnomad_path
    tb = pysam.TabixFile(db)
       
    with open(input_file, 'r') as hin:
        with open(output_file,'w') as hout:
            for line in hin: # one SJ
                line = line.rstrip('\n')
                junc_record = line
                F = line.split('\t')
                if F[23] != "True": continue

                c = (F[0].replace('chr', '')) 
   
                if c != "Y": 
                # if c in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]:
                    if F[13] == "-":
                        out_record = line + "\t-\t-\t-\t-\t0.0\n" 
                        hout.write(out_record)
                    else:
                        out_record = ""
                        if genome_id == "hg19":
                            chr = str(c)
                        else:
                            chr = "chr" + str(c)

                        rows = tb.fetch(chr, int(F[16]) - 1, int(F[16]))

                        cur_AF = 0.0
                        cur_allele = "-"
                        if rows is not None:
                            for row in rows:
                                srow=str(row)
                                record = srow.split('\t')
                                 
                                if F[17] == record[3] and F[18] == record[4]:
                                    allele = record[3]+">"+record[4]
                                    infos = record[7].split(';')
                                    for info in infos:
                                        if info.startswith("AF="):
                                            cur_AF = float(info.replace("AF=", ''))
                                            cur_allele = allele

                                    break                                        
           
                        out_record = junc_record + "\t" + cur_allele + "\t" + str(cur_AF) +"\n" 
                        hout.write(out_record)
    
                                    
                elif c in ["Y"]:
                    out_record = junc_record + "\tna\t0\n" 
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
