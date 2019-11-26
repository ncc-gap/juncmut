#! /usr/bin/env python

import os
import subprocess
import pandas as pd
import pysam
import re

def run_juncm(args):
    import shutil
    from junc_utils import utils
    if not os.path.exists('alterativeSJ_fil_annot'):
        os.mkdir('alterativeSJ_fil_annot')
    if not os.path.exists('alterativeSJ_assadjfreq'):
        os.mkdir('alterativeSJ_assadjfreq')
    if not os.path.exists('alterativeSJ_tabixAI'):
        os.mkdir('alterativeSJ_tabixAI')
    if not os.path.exists('alterativeSJ_cmut'):
        os.mkdir('alterativeSJ_cmut')
    
    basename = os.path.basename(args.input_file)  #./junction/*.SJ.out.tab    
    pr=basename.split('.', 4)[0]    
    print(pr)    
    file1 = './alterativeSJ_fil_annot/%s.SJ.fil.txt' %pr
    file2 = './alterativeSJ_fil_annot/%s.SJ.fil.annot.txt' %pr
    
    if not args.control_file:
        shutil.copy(args.input_file, file1)

    else:
        f = args.input_file
        cont_list=args.control_file
        n=1
        for cont in cont_list:
            out = 'tmp_out'+str(n)+'.txt'
            #junc_utils filter --pooled_control_file *.bed.gz input(./junction/A427.SJ.out.tab) output
            utils.proc_star_junction(f, out, cont, 3, 10, True, False)
            f = 'tmp_in'+str(n)+'.txt'
            shutil.copy(out, f)
        n=n+1
    
    shutil.copy(out, file1)
    
    #junc_utils annotate ${out2} ${out3} --genome_id hg19
    #annotate.annot_junction(args.junc_file, args.output_path, args.junction_margin, args.exon_margin, args.genome_id, is_grc) 
    #def annotate_main(args): 
    #from junc_utils import annotate
    #from annot_utils.utils import grc_check
    
    #if args.grc == True:
    #    logger.warning("--grc argument is deprecated and ignored.")

    #is_grc = grc_check(args.junc_file, [0])
 
    #annotate.annot_junction(args.junc_file, args.output_path, args.junction_margin, args.exon_margin, args.genome_id, is_grc)
    #annotate.annot_junction(output_file1, output_file2, junction_margin=3, exon_margin=30, genome_id="hg19") #, is_grc) #is_grc=False
    annotate_commands = ["junc_utils", "annotate", file1, file2, "--genome_id", args.genome_id]
    subprocess.call(annotate_commands)
    
    print("Filtering control and annotation were completed.")
###ass_adj_freq###    
#input ./alterativeSJ_filtered_annot/sample.SJ.fil.annot.txt
#output ./alterativeSJ_adjustment/sample.SJ.fil.annot.ass.adj.freq.txt
#Function
#extract alterative SS (ass.
#ajustment of intron start end (adj.
#calsulate frequency

#def asj_ass_adj_bed(input_file, output_file, pr):
    file3='./alterativeSJ_assadjfreq/%s.SJ.fil.annot.assadj.txt' %pr
    with open(file2, 'r') as in1:
        with open(file3, 'w') as out1:

            for line in in1:
                F = line.rstrip('\n').split('\t')
                r=F[6]
                if F[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
                    if "Alternative" in F[9] or "alternative" in F[9]:
                        if "e" in F[13]: #strand=+ 3'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+'\t'+str(si)+'\t'+str(ei)+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "s" in F[13]: #strand=- 5'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                        elif "s" in F[17]: #strand=+ 5'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "e" in F[17]: #strand=- 3'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(pr)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')

#def asj_freq(input_file, output_file, pr):
   
    jfile= './junction/%s.SJ.out.tab' %(pr)
    jf = pd.read_csv(jfile, sep='\t', header=None, index_col=None,names=('CHR','START','END','A','B','C','reads','D','E'),
                 dtype = {'CHR':'object', 'START':'object', 'END':'object','A':'object','B':'object','C':'object','reads':'int','D':'object','E':'object'})
    file4 ='./alterativeSJ_assadjfreq/%s.SJ.fil.annot.assadjfreq.txt' %pr 
    with open(file3, 'r') as df1:
        with open(file4, 'w') as out1:
            for line in df1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                if '5' in F[6] and '+' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[4], F[2])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['END'].isin(pl))]  
                    #total = rec.groupby('END').sum()['reads'] <--XX
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
                elif '5' in F[6] and '-' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[3], F[1])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['START'].isin(pl))]  
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
                elif '3' in F[6] and '+' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[3], F[1])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['START'].isin(pl))]  
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
                elif '3' in F[6] and '-' in F[7]:
                    c = (F[0],)
                    cl = list(c)
                    p = (F[4], F[2])
                    pl = list(p)
                    rec = jf[(jf['CHR'].isin(cl) ) & (jf['END'].isin(pl))]  
                    total = int(rec.sum()['reads'])
                    freq = round(int(F[8])/total, 3)
                    out1.write(ln+'\t'+str(total)+'\t'+str(freq)+'\n')
	
    print("Extraction of SS was completed.")
    print("Frequency culculatetion was completed.")
#def asj_tabix_hg19(input_file, output_file, pr):

    tbx = pysam.TabixFile("./reference/whole_genome_filtered_spliceai_scores.vcf.gz")
    file5 = './alterativeSJ_tabixAI/%s.SJ.fil.annot.assadjfreq.AI.txt' %(pr)
    with open(file4, 'r') as in1:
        with open(file5, 'w') as out1:
            for line in in1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                
                if F[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
                    if "3" in F[6] and "+" in F[7]:
                        c=F[0]
                        s=int(F[2])-5
                        e=int(F[2])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            #strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if inte == aip+ag-1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
                    elif "3" in F[6] and "-" in F[7]:
                        c=F[0]
                        s=int(F[1])-5
                        e=int(F[1])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            #strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if ints == aip+ag+1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
                    elif "5" in F[6] and "-" in F[7]:
                        c=F[0]
                        s=int(F[2])-5
                        e=int(F[2])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            #strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if inte == aip+dg-1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
                    elif "5" in F[6] and "+" in F[7]:
                        c=F[0]
                        s=int(F[1])-5
                        e=int(F[1])+5
                        rows = tbx.fetch(c, s, e)
                        num=0
                        for row in rows:
                            line=ln+'\t'+str(row)
                            G = line.split('\t')
                            cl = G[4] #class
                            ints = int(G[1]) #SAVnet intron start
                            inte = int(G[2])
                            aip = int(G[12]) #Splice AI postion
                            inf = G[18]
                            inf2 = re.split('[=;]', inf)
                            #strand = G[7]
                            ag = int(inf2[17])
                            dg = int(inf2[21])
                            del G[-1]
                            line= '\t'.join(G)
                            G_list=(inf2[1],inf2[3],inf2[5],inf2[7],inf2[9],inf2[11],inf2[13],inf2[15],inf2[17],inf2[19],inf2[21],inf2[23])
                            line2= '\t'.join(G_list)
                            if ints == aip+dg+1:
                                num+=1
                                out1.write(line+'\t'+line2+'\tc'+str(num)+'\n')
                            else:
                                out1.write(line+'\t'+line2+'\t-\n')
 
    print("Annotation of spliceAI was completed.")
#def asj_cmut #extraction of_canonical mutation
    
    def complement_dna(string):
        comp=''
        for char in string:
            if   char == 'A': comp += 'T'
            elif char == 'T': comp += 'A'
            elif char == 'G': comp += 'C'
            elif char == 'C': comp += 'G'
            else:             comp += char
        return comp[::-1]

    
    file6 = './alterativeSJ_cmut/%s_fil.annot.assadjfreq.AI.cmut.txt' %(pr)
    he = 'chr\tstart\tend\tstart_ori\tend_ori\tsample\tclass\tstrand\treads\ttotal\tfreq\tCHR\tPOS\tID\tREF\tMUT\tSCORE1\tSCORE2\tSYMBOL\tSTRAND\tTYPE\tDIST\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\tSAVnet_SpliceAI\tMOTIF\tDP\tDS\n'
    with open(file6, 'w') as out1:
        out1.write(he)
    
        with open(file5, 'r') as in1:   
        
            for line in in1:
                F = line.rstrip('\n').split('\t')
                ln = line.rstrip('\n')
                if int(F[8])>= args.read_num_thres and float(F[10])>= args.freq_thres and "c" in F[30]:
        #3'SS strand=+
                    if "3" in F[6] and "+" in F[7]:
                        a =str(F[22])
                        DP =str(F[26])
                        if F[26] == "1" and F[15] == "G":
                            m= "*"+F[14]+"|_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[26] == "2" and F[15] == "A":
                            m= F[14]+"*|_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #3'SS strand=- 
                    elif "3" in F[6] and "-" in F[7]:
                        rb=complement_dna(F[14])
                        mb=complement_dna(F[15])
                        a =str(F[22])
                        DP =str(F[26])
                        if F[26] == "-1" and F[15] == "C":
                            m= "*"+rb+"|_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[26] == "-2" and F[15] == "T":
                            m= rb+"*|_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #5'SS strand=+
                    elif "5" in F[6] and "+" in F[7]:
                        a =str(F[24])
                        DP =str(F[28])
                        if F[28] == "-1" and F[15] == "G":
                            m= "|"+F[14]+"*_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[28] == "-2" and F[15] == "T":
                            m= "|*"+F[14]+"_"+F[15]
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
        #5'SS strand=-
                    elif "5" in F[6] and "-" in F[7]:
                        rb=complement_dna(F[14])
                        mb=complement_dna(F[15])
                        a =str(F[24])
                        DP =str(F[28])
                        if F[28] == "1" and F[15] == "C":
                            m= "|"+rb+"*_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')
                        elif F[28] == "2" and F[15] == "A":
                            m= "|*"+rb+"_"+mb
                            out1.write(ln+'\t'+m+'\t'+str(DP)+'\t'+str(a)+'\n')

    print("The prediction of canonical mutations was completed.")
