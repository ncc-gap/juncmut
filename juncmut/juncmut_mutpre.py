#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 2019

@author: genome

#def asj_mutcall_hg19(input_pr, folder):

"""

def juncmut_mutpre(input_file, output_file, reference):

    import pysam
        
    ref = pysam.FastaFile(reference) # 0-coordinate
    
    trans = str.maketrans('ATGCatgc', 'TACGTAGC')
    
    dict5p = {0: 'A', 1: 'G', 2: 'G', 3: 'T', 4: 'AG', 5: 'A', 6: 'G', 7: 'T'}
    dict3p = {0: 'CT', 1: 'CT', 2: 'ATGC', 3: 'CT', 4: 'A', 5: 'G', 6: 'AG'}
    dict5m = {0: 'A', 1: 'C', 2: 'T', 3: 'CT', 4: 'A', 5: 'C', 6: 'C', 7: 'T'}
    dict3m = {0: 'CT', 1: 'C', 2: 'T', 3: 'AG', 4: 'TCGA', 5: 'AG', 6: 'AG'}
    
    
    with open(input_file, 'r') as in1:
        with open(output_file, 'w') as out1:
            for line in in1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                
                if "5" in F[6] and "+" in F[7]: #o-->
                    c=F[0]
                    #c = (c.replace('chr', ''))
                    s=int(F[1])-1-2 #region -2|0~5  
                    e=int(F[1])+5
                    getbases = ref.fetch(c,s,e)
                    bases = getbases.upper()

                    mlist = []
                    mposi = []
                    noncamposi = []
                    for i in range(0,len(bases)): 
                        if bases[i] in dict5p[i]:
                            mlist.append(".")
                            mposi.append(str(int(F[1])+int(i-2))+':'+bases[i]+':'+dict5p[i])
                        else:
                            mlist.append(dict5p[i])
                            mposi.append(str(int(F[1])+int(i-2))+':'+bases[i]+':'+dict5p[i])
                            noncamposi.append(str(int(F[1])+int(i-2))+':'+bases[i]+':'+dict5p[i])
                    mlist[4] = (mlist[4].replace('AG', 'R'))
                    
                    if mlist[2] == "G" or mlist[3] == "T":
                        
                        if mlist[2] == "G" and mlist[3] == "T":
                            mposi_rec = mposi[2] + "," + mposi[3]
                            allele = '|' + bases[2] + bases[3] + '>GT*d'
                        elif mlist[2] == "G" and mlist[3] == ".":
                            mposi_rec = mposi[2]
                            allele = '|' + bases[2] + bases[3] + '>GT'
                        elif mlist[2] == "." and mlist[3] == "T":
                            mposi_rec = mposi[3]
                            allele = '|' + bases[2] + bases[3] + '>GT'
                    elif len(set(mlist)) == 1:
                        allele = '-'
                        mposi_rec = '-'
                    else:
                        allele = 'non-canonical'
                        mposi_rec = ','.join(map(str,noncamposi))                        
                        
                    rec = ln + "\t" + bases + "\t" + ''.join(mlist) + "\t" +  str(mposi_rec) + "\t"+ allele + "\n"
                    out1.write(rec)             
                
                elif "5" in F[6] and "-" in F[7]: #<--o
                    c=F[0]
                    #c = (c.replace('chr', ''))
                    s=int(F[2])-1-5
                    e=int(F[2])+2
                    getbases = ref.fetch(c,s,e)
                    bases = getbases.upper()
                    #rev = ''.join(reversed(bases))
                    #comp = rev.translate(trans)
                    
                    mlist = []
                    mposi = []
                    noncamposi = []
                    for i in range(0,len(bases)):
                        if bases[i] in dict5m[i]:
                            mlist.append(".")
                            mposi.append(str(int(F[2])+int(i-5))+':'+bases[i]+':'+dict5m[i])
                        else:
                            mlist.append(dict5m[i])
                            mposi.append(str(int(F[2])+int(i-5))+':'+bases[i]+':'+dict5m[i])
                            noncamposi.append(str(int(F[2])+int(i-5))+':'+bases[i]+':'+dict5m[i])
                    mlist[3] = (mlist[3].replace('CT', 'Y'))
                    
                    if mlist[4] == "A" or mlist[5] == "C":

                        if mlist[4] == "A" and mlist[5] == "C":
                            mposi_rec = mposi[4] + "," + mposi[5]
                            allele = '|' + bases[5].translate(trans) + bases[4].translate(trans) + '>GT*d'
                        elif mlist[4] == "A" and mlist[5] == ".":
                            mposi_rec = mposi[4]
                            allele = '|' + bases[5].translate(trans) + bases[4].translate(trans) + '>GT'
                        elif mlist[4] == "." and mlist[5] == "C":
                            mposi_rec = mposi[5] 
                            allele = '|' + bases[5].translate(trans) + bases[4].translate(trans) + '>GT'
                    elif len(set(mlist)) == 1:
                        allele = '-'
                        mposi_rec = '-'
                    else:
                        allele = 'non-canonical' 
                        mposi_rec = ','.join(map(str,noncamposi))
                    rec = ln + "\t" + bases + "\t" + ''.join(mlist) + "\t" +  str(mposi_rec) + "\t"+ allele + "\n"
                    out1.write(rec) 
                    
                elif "3" in F[6] and "+" in F[7]: #-->o
                    c=F[0]
                    #c = (c.replace('chr', ''))
                    s=int(F[2])-1-5 #region -5|1  
                    e=int(F[2])+1
                    getbases = ref.fetch(c,s,e)
                    bases = getbases.upper()
                    mlist = []
                    mposi = []
                    noncamposi = []
                    for i in range(0,len(bases)):
                        if bases[i] in dict3p[i]:
                            mlist.append(".")
                            mposi.append(str(int(F[2])+int(i-5))+':'+bases[i]+':'+dict3p[i])
                        else:
                            mlist.append(dict3p[i])
                            mposi.append(str(int(F[2])+int(i-5))+':'+bases[i]+':'+dict3p[i])
                            noncamposi.append(str(int(F[2])+int(i-5))+':'+bases[i]+':'+dict3p[i])
                    mlist[0] = (mlist[0].replace('CT', 'Y'))
                    mlist[1] = (mlist[1].replace('CT', 'Y'))
                    mlist[3] = (mlist[3].replace('CT', 'Y'))
                    mlist[6] = (mlist[6].replace('AG', 'R'))

                    if mlist[4] == "A" or mlist[5] == "G":
                        if mlist[4] == "A" and mlist[5] == "G":
                            mposi_rec = mposi[4] + "," + mposi[5]
                            allele =  bases[4] + bases[5] + '|>AG*d'
                        elif mlist[4] == "A" and mlist[5] == ".":
                            mposi_rec = mposi[4]
                            allele =  bases[4] + bases[5] + '|>AG'
                        elif mlist[4] == "." and mlist[5] == "G":
                            mposi_rec = mposi[5]
                            allele =  bases[4] + bases[5] + '|>AG'
                        
                    elif len(set(mlist)) == 1:
                        allele = '-'
                        mposi_rec = '-'
                    else:
                        allele = 'non-canonical' 
                        mposi_rec = ','.join(map(str,noncamposi))  
                    rec = ln + "\t" + bases + "\t" + ''.join(mlist) + "\t" +  str(mposi_rec) + "\t"+ allele + "\n"
                    out1.write(rec)              
                
                elif "3" in F[6] and "-" in F[7]: #o<--
                    c=F[0]
                    #c = (c.replace('chr', ''))
                    s=int(F[1])-1-1
                    e=int(F[1])+5
                    getbases = ref.fetch(c,s,e)
                    bases = getbases.upper()
                    #rev = ''.join(reversed(bases))
                    #comp = rev.translate(trans)
                    mlist = []
                    mposi = []
                    noncamposi = []
                    for i in range(0,len(bases)):
                        if bases[i] in dict3m[i]:
                            mlist.append(".")
                            mposi.append(str(int(F[1])+int(i-1))+':'+bases[i]+':'+dict3m[i])
                        else:
                            mlist.append(dict3m[i])
                            mposi.append(str(int(F[1])+int(i-1))+':'+bases[i]+':'+dict3m[i])
                            noncamposi.append(str(int(F[1])+int(i-1))+':'+bases[i]+':'+dict3m[i])
                    mlist[0] = (mlist[0].replace('CT', 'Y'))
                    mlist[3] = (mlist[3].replace('AG', 'R'))
                    mlist[5] = (mlist[5].replace('AG', 'R'))
                    mlist[6] = (mlist[6].replace('AG', 'R'))
                    

                    if mlist[1] == "C" or mlist[2] == "T":
                        if mlist[1] == "C" and mlist[2] == "T":
                            mposi_rec = mposi[1] + "," + mposi[2]
                            allele =  bases[2].translate(trans) + bases[1].translate(trans) + '|>AG*d'
                        elif mlist[1] == "C" and mlist[2] == ".":
                            mposi_rec = mposi[1]
                            allele =  bases[2].translate(trans) + bases[1].translate(trans) + '|>AG'
                        elif mlist[1] == "." and mlist[2] == "T":
                            mposi_rec = mposi[2]
                            allele =  bases[2].translate(trans) + bases[1].translate(trans) + '|>AG'                        
                        
                    elif len(set(mlist)) == 1:
                        allele = '-'
                        mposi_rec = '-'
                    else:
                        allele = 'non-canonical' 
                        mposi_rec = ','.join(map(str,noncamposi)) 
                    rec = ln + "\t" + bases + "\t" + ''.join(mlist) + "\t" +  str(mposi_rec) + "\t"+ allele + "\n"
                    out1.write(rec)  

if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input_file", metavar = "input_file", default = None, type = str,
                            help = "prefix") 
        
    parser.add_argument("output_file", metavar = "output_file", default = "my_samples", type = str,
                            help = "folder name of input files") 
    
    parser.add_argument("reference", metavar = "reference", default = "my_samples", type = str,
                            help = "/full/path/to/reference")     
        
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    reference = args.reference
    
    juncmut_mutpre(input_file, output_file, reference)
