#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 09:25:52 2020

@author: genome

USAGE
python juncmut_postProcess.py input_file GTEx_control_file gene_file cancer_gene_file
"""

def postProcess(input_file):

    import pathlib

    in_path = pathlib.Path(input_file)
    output_file = in_path.stem+'.motif.txt'
    
    
    hout = open(output_file, 'w') 
    header = ["Sample", "Chr", "SJ_Start", "SJ_End", "SJ_Type", "SJ_Strand", "SJ_Read_Count", "SJ_Depth", "SJ_Freq",
                  "Ref_Motif", "Possivle_Alt_Motif", "Is_Canonical", "SJ_Overlap_Count", 
                  "Mut_Pos", "Mut_Ref", "Mut_Alt", "Mut_Count", "Mut_Depth", "Mut_Freq", "AF_Gnomad",
                  "Realign_No_SJ_Neg", "Realign_No_SJ_Pos", "Realign_Target_SJ_Neg", "Reaglin_Target_SJ_Pos",
                  "Realign_Normal_SJ_Neg", "Realign_Normal_SJ_Pos","SUM(Realign_No_SJ_Pos+Reaglin_Target_SJ_Pos)",
                  "RATIO(Reaglin_No_SJ_Pos/Realign_No_SJ)","SJ_intronic_type","Mut_Motif_Pos"]
    
    print('\t'.join(header), file = hout)
    hout.close()
    
    
    with open(input_file, 'r') as in1, open(output_file, 'a') as hout:
        next(in1)
        for line in in1:
            F = line.rstrip('\n').split('\t')
            ln = line.rstrip('\n')
            
            sum_vx = int(F[21])+int(F[23])
            sum_uv = int(F[20])+int(F[21])
            if sum_uv == 0:
                ratio_uv = 'nan'
            else:
                ratio_uv = int(F[21])/(int(F[20])+int(F[21]))
            
            sj_type = str(F[4])
            if 'Intronic' in sj_type:
                rec_sj_type = 'Intronic'
            else:
                rec_sj_type = 'Neighboring_exonic/Exonic'
                
            #o-->
            if "5" in F[4] and "+" in F[5]: 
                dis_mut_SJstart = int(F[13])-int(F[2])
                if dis_mut_SJstart <0 and dis_mut_SJstart >-3:
                    mutation = 'exon'
                elif dis_mut_SJstart <6 and dis_mut_SJstart >=0:
                    mutation = 'intron'
                else:
                    mutation = 'outside of motif'
            #-->o
            elif "3" in F[4] and "+" in F[5]: 
                dis_mut_SJend = int(F[13])-int(F[3])
                if dis_mut_SJend ==1:
                    mutation = 'exon'
                elif dis_mut_SJend <1 and dis_mut_SJend >-6:
                    mutation = 'intron'
                else:
                    mutation = 'outside of motif'
            #o<--       
            elif "3" in F[4] and "-" in F[5]: 
                dis_mut_SJstart = int(F[13])-int(F[2])
                if dis_mut_SJstart ==-1:
                    mutation = 'exon'
                elif dis_mut_SJstart <6 and dis_mut_SJstart >=0:
                    mutation = 'intron'
                else:
                    mutation = 'outside of motif'
            #<--o
            elif "5" in F[4] and "-" in F[5]: 
                dis_mut_SJend = int(F[13])-int(F[3])
                if dis_mut_SJend >0 and dis_mut_SJend <3:
                    mutation = 'exon'
                elif dis_mut_SJend <1 and dis_mut_SJend >-6:
                    mutation = 'intron'
                else:
                    mutation = 'outside of motif'
            rec = ln +'\t'+ str(sum_vx) +'\t'+ str(ratio_uv) +'\t'+ rec_sj_type +'\t'+ mutation +'\n'
            hout.write(rec)

def add_GTEx_control(input_file, cont):
    import pysam
    import pathlib
    
    db = cont      
    tb = pysam.TabixFile(db)        

    in_path = pathlib.Path(input_file)
    output_file = in_path.stem+'.gtex.txt'
    
    hout = open(output_file, 'w') 
    header = ["Sample", "Chr", "SJ_Start", "SJ_End", "SJ_Type", "SJ_Strand", "SJ_Read_Count", "SJ_Depth", "SJ_Freq",
              "Ref_Motif", "Possivle_Alt_Motif", "Is_Canonical", "SJ_Overlap_Count", 
              "Mut_Pos", "Mut_Ref", "Mut_Alt", "Mut_Count", "Mut_Depth", "Mut_Freq", "AF_Gnomad",
              "Realign_No_SJ_Neg", "Realign_No_SJ_Pos", "Realign_Target_SJ_Neg", "Reaglin_Target_SJ_Pos",
              "Realign_Normal_SJ_Neg", "Realign_Normal_SJ_Pos","SUM(Realign_No_SJ_Pos,Reaglin_Target_SJ_Pos)",
              "SUM(Realign_No_SJ_Pos+Reaglin_Target_SJ_Pos)","SJ_Int_exonic_type","Mut_Motif_Pos","GTEx_info","GTEx_count"]

    print('\t'.join(header), file = hout)

    with open(input_file) as hin:
        next(hin)
        for line in hin:
            F = line.rstrip('\n').split('\t')
            ln = line.rstrip('\n')
            c = F[1]
            s = int(F[2])
            e = int(F[3])
        
            rows =tb.fetch(c, s, e)
            
            sj_list = []
            count_list =[]
            if rows is not None:
                for row in rows:
                    record = str(row).split('\t')
                    if abs(s-int(record[1]))<5 and abs(e-int(record[2]))<5:
                        sj_list.append(record[3])
                        count = len(record[3].split(','))
                        count_list.append(str(count))
            
            if sj_list:
                sj_out = ';'.join(sj_list)
                count_out = ';'.join(count_list)
            else:
                sj_out = '-'
                count_out = '-'
            record_out = ln +"\t"+ str(sj_out) +"\t"+ str(count_out) + "\n" 
            hout.write(record_out)
            
    hout.close()

def add_gene(input_file, gene_file):
    
    import pysam
    import pathlib
    
    gene_tb = pysam.TabixFile(gene_file)

    in_path = pathlib.Path(input_file)
    output_file = in_path.stem+'.gene.txt'
    
    hout = open(output_file, 'w') 
    
    with open(input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        
        print('\t'.join([header[0], header[1], header[2], header[3]] + ["Gene"] + header[4:]), file = hout)
    
        for line in hin:
            F = line.rstrip('\n').split('\t')
    
            tchr = 'chr' + F[1]
            tpos_s = F[2]
            tpos_e = F[3]
    
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch(tchr, int(tpos_s)-5, int(tpos_e) + 5)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorFlag = 1
    
            gene = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene.append(record[3])
    
            if len(gene) == 0: gene.append("---")
            gene = ';'.join(list(set(gene)))
    
            print('\t'.join([F[0], F[1],F[2],F[3]] + [gene] + F[4:]), file = hout)
            
    hout.close()


def add_cancer_gene(input_file, cgene_file):
    
    import pathlib

    in_path = pathlib.Path(input_file)
    output_file = in_path.stem+'.cgene.txt'
    
    cancer_dict = {}
    with open(cgene_file, 'r') as db:
        for data in db:
            D = data.rstrip('\n').split('\t')
            key = D[1]
            val = '\t'.join([D[1],D[3],D[4],D[5],D[6]])
            
            cancer_dict[key] = val
   
  
    
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        
        header = hin.readline().rstrip('\n').split('\t')
        hout.write('\t'.join(header[0:]) +"\tCancer_gene\tElement_type\tCategory\tMoF\tTissue\n")      
    
        for line in hin:

            F = line.rstrip('\n').split('\t')
    
            q_gene = F[4]  
            v_gene = cancer_dict.get(q_gene, '-\t-\t-\t-\t-')                    
            rec = '\t'.join(F[0:]+ [v_gene]) +'\n'
            hout.write(rec)
            


if __name__== "__main__":
    import argparse
    import pathlib
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input_file", metavar = "input_file", default = None, type = str,
                            help = "Path to input file") 
    
    parser.add_argument("cont", metavar = "cont", default = "gtex_sj2r2s.bed.gz", type = str,
                            help = "Path to control file") 
        
    parser.add_argument("gene_file", metavar = "gene_file", default = "wgEncodeGencodeBasicV33lift37_hg19.bed.gz", type = str,
                            help = "Path to gene file") 
    
    parser.add_argument("cgene_file", metavar = "cgene_file", default = "TableS1_compendium_mutational_drivers.txt", type = str,
                            help = "Path to cancer gene file") 
        
    args = parser.parse_args()
    
    input_file = args.input_file
    cont = args.cont
    gene_file = args.gene_file
    cgene_file = args.cgene_file
    
    #input_file = 'suzuki_cellline.juncmut.result.val.txt'

    in_path = pathlib.Path(input_file)
    input_file2 = in_path.stem+'.motif.txt'
    input_file3 = in_path.stem+'.motif.gtex.txt'
    input_file4 = in_path.stem+'.motif.gtex.gene.txt'
    
    postProcess(input_file)
    
    add_GTEx_control(input_file2, cont)
    
    add_gene(input_file3, gene_file)
    
    add_cancer_gene(input_file4, cgene_file)
         
