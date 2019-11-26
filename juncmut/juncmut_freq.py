#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Naoko Iida
Date:10152019

"""

def juncmut_freq(pr, folder, read_num_thres, freq_thres):
    import pandas as pd
    
    file3='./data/%s/alterativeSJ_assadjfreq/%s.SJ.fil.annot.assadj.txt' %(folder,pr)
    jfile= './data/%s/%s.SJ.out.tab' %(folder,pr)
    file4 ='./data/%s/alterativeSJ_assadjfreq/%s.SJ.fil.annot.assadjunifreq.txt' %(folder,pr) 
    file44 ='./data/%s/alterativeSJ_assadjfreq/%s.SJ.fil.annot.assadjunifreqT.txt' %(folder,pr) 
      
    data = pd.read_csv(file3, sep='\t', header=None,  dtype={0:'str',1:'int',2:'int',3:'int',4:'int'})
    data.columns = ['chr','s','e','s_ori','e_ori', 'sample', 'class','strand', 'reads']
    data['junc'] = data[['chr','s','e']].apply(lambda x: '{}:{}:{}'.format(x[0],x[1],x[2]), axis=1)
    group_junc = data.groupby(['junc'])
    agg_junc = group_junc.agg({"chr": "max", "s": "unique", "e": "unique", "s_ori": "unique","e_ori": "unique", 'sample':'unique', 'class':'unique', 'strand':'unique', 'reads':'sum'}) #"reads": "max", 
    list_junc = agg_junc.sort_values(by=["junc"], ascending=False) #data.frame
                         
    with open(jfile) as d1:
        data1 = pd.read_csv(d1, delimiter='\t',usecols=[0,1,2,6], header=None, dtype={0:'str',1:'int',2:'int',6:'int'}) 
        data1.columns = ['chr', 'start','end', 'reads']
        s_reads= data1.groupby(['chr','start'])['reads'].sum()
        s_dict = s_reads.to_dict()
        e_reads= data1.groupby(['chr','end'])['reads'].sum()
        e_dict = e_reads.to_dict()
        
    with open(file4, 'w') as out1:
        for row in list_junc.itertuples():
    
            junc = str(row[0])
            c = str(row[1])
            
            s = set(row[2])
            s_ori = set(row[4])
            start= s | s_ori
            start_l = list(start)
            start_ls = sorted(start_l)
            
            e = set(row[3])
            e_ori = set(row[5])
            end= e | e_ori
            end_l = list(end)
            end_ls = sorted(end_l)   
            
            total = 0
                #strand=+ 5'SS end-side
            if "5" in str(row[7]) and "+" in str(row[8]):  
                for i in range(0,len(end_ls)):
                    
                    v=e_dict.get((str(c), int(end_ls[i])), '0')
                    total = total + int(v) 
                freq = int(row[9])/total
                rec = c + "\t" + str(''.join(map(str, s))) + "\t" + str(''.join(map(str, e))) + "\t"  + str(';'.join(map(str, start_ls))) + "\t" +  str(';'.join(map(str, end_ls)))  + "\t" + \
                str(''.join(map(str, row[6]))) + "\t" +  str(''.join(map(str, row[7]))) + "\t" + str(''.join(map(str, row[8]))) + "\t" + str(row[9]) + "\t" + str(total)+ "\t" + str(freq) + '\n' #depth
                out1.write(rec)
               #strand=- 3'SS end-side                                   
            elif "3" in str(row[7]) and "-" in str(row[8]):
                for i in range(0,len(end_ls)):
                    v=e_dict.get((str(c), int(end_ls[i])), '0')
                    total = total + int(v) 
                freq = int(row[9])/total
                rec = c + "\t" + str(''.join(map(str, s))) + "\t" + str(''.join(map(str, e))) + "\t"  + str(';'.join(map(str, start_ls))) + "\t" +  str(';'.join(map(str, end_ls)))  + "\t" + \
                str(''.join(map(str, row[6]))) + "\t" +  str(''.join(map(str, row[7]))) + "\t" + str(''.join(map(str, row[8]))) + "\t" + str(row[9]) + "\t" + str(total)+ "\t" + str(freq) + '\n' #depth
                out1.write(rec)                
                
                #strand=- 5'SS start-side         
            elif "5" in str(row[7]) and "-" in str(row[8]):
                for i in range(0,len(start_ls)):
                    v=s_dict.get((str(c), int(start_ls[i])), '0')
                    total = total + int(v) 
                freq = int(row[9])/total
                rec = c + "\t" + str(''.join(map(str, s))) + "\t" + str(''.join(map(str, e))) + "\t"  + str(';'.join(map(str, start_ls))) + "\t" +  str(';'.join(map(str, end_ls)))  + "\t" + \
                str(''.join(map(str, row[6]))) + "\t" +  str(''.join(map(str, row[7]))) + "\t" + str(''.join(map(str, row[8]))) + "\t" + str(row[9]) + "\t" + str(total)+ "\t" + str(freq) + '\n' #depth
                out1.write(rec)                
                        
             #strand=+ 3'SS start-side        
            elif "3" in str(row[7]) and "+" in str(row[8]):
                for i in range(0,len(start_ls)):
                    v=s_dict.get((str(c), int(start_ls[i])), '0')
                    total = total + int(v) 
                freq = int(row[9])/total
                rec = c + "\t" + str(''.join(map(str, s))) + "\t" + str(''.join(map(str, e))) + "\t"  + str(';'.join(map(str, start_ls))) + "\t" +  str(';'.join(map(str, end_ls)))  + "\t" + \
                str(''.join(map(str, row[6]))) + "\t" +  str(''.join(map(str, row[7]))) + "\t" + str(''.join(map(str, row[8]))) + "\t" + str(row[9]) + "\t" + str(total)+ "\t" + str(freq) + '\n' #depth
                out1.write(rec)  
                
    with open(file4, 'r') as in1:
            with open(file44, 'w') as out2:
                for line in in1:
                    F = line.rstrip('\n').split('\t')
                    if int(float(F[8])) >= read_num_thres and float(F[10]) >= freq_thres:
                    #if int(F[8]) >= 3 and float(F[10]) >= 0.05:
                        out2.write(line)
                    
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input", metavar = "prefix", default = None, type = str,
                            help = "Path to input file") 
        
    parser.add_argument("folder", metavar = "group", default = "samples", type = str,
                            help = "Path to input file") 
    
    parser.add_argument("--read_num_thres", type = int, default = 3,
                            help = "Remove splicing junctions whose supporting numbers are below this value (default: %(default)s)")
        
    parser.add_argument("--freq_thres", type = float, default = 0.05,
                            help = "Remove splicing junctions whose frequency is below this value (default: %(default)s)")
    
    args = parser.parse_args()
    
    pr = args.input
    folder = args.folder
    read_num_thres = args.read_num_thres
    freq_thres = args.freq_thres
    
    juncmut_freq(pr, folder, read_num_thres, freq_thres)      
