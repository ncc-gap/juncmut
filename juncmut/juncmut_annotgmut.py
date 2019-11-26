#! /usr/bin/env python
"""
Naoko Iida
python ./juncagg_codes/juncmut_annotgmut.py folder_name input_file (folder and file under juncmut_agg folder)
"""
#import pandas as pd
import subprocess
import pandas as pd
import os
import pathlib
import csv


def tidy_bases(bases, qualities):  #remove indel and edeage
    
    import re
    #bases = ",,.,^!A$a^-A^!a"
    #qualities = "EEDE@D:C"
    proc1 = ""
    proc2 = ""
    
    #del end position in a read.   $ is the last position of read, following the base. ^ is the start position of read, following the quality and the base.           
    while len(bases) > 0:
        match = re.search(r'[$\^]', bases)
        if match is None:
            proc1 = proc1 + bases
            bases = ""
            proc2 = proc2 + qualities
            qualities = ""
        elif match == "^":
            pos = match.start()
            proc1 = proc1 + bases[0:pos]
            bases = bases[(pos + 3):len(bases)]
            proc2 = proc2 + qualities[0:pos]
            qualities = qualities[(pos + 1): len(qualities)] #delete 1 char
        else:
            pos = match.start()
            proc1 = proc1 + bases[0:pos]
            bases = bases[(pos + 2):len(bases)]
            proc2 = proc2 + qualities[0:pos]
            qualities = qualities[(pos + 1): len(qualities)] #delete 1 char
    #del indel. +/- +[0-9]+[bases]. * means the deletion(?).
    bases = proc1
    proc1 = ""
    qualities = proc2
    proc2 = ""
    while len(bases)>0:
        match = re.search(r'[\+\-\*]', bases)
        if match is None:
            proc1 = proc1 + bases
            bases = ""
            proc2 = proc2 + qualities
            qualities = ""
        elif match.group() == "*":
            pos = match.start()
            proc1 = proc1 + bases[0:pos]
            bases = bases[pos+1: len(bases)]
            proc2 = proc2 + qualities[0:pos]
            qualities = qualities[(pos + 1): len(qualities)]     
        else:
            pos = match.start()
            del_num = bases[match.start()+1]
            proc1 = proc1 + bases[0:pos]
            bases = bases[pos+int(del_num)+2: len(bases)]
            proc2 = proc2 + qualities[0:pos]
            qualities = qualities[(pos + 1): len(qualities)] #delete 1 char
    #del no-base position in a read. 
    bases = proc1
    proc1 = ""
    qualities = proc2
    proc2 = ""
    while len(bases) > 0:
        match = re.search(r'[<>]', bases)
        if match is None:
            proc1 = proc1 + bases
            bases = ""
            proc2 = proc2 + qualities
            qualities = ""
        else:
            pos = match.start()
            proc1 = proc1 + bases[0:pos]
            bases = bases[(pos + 1):len(bases)]
            proc2 = proc2 + qualities[0:pos]
            qualities = qualities[(pos + 1): len(qualities)]

    #del low quality base
    bases = proc1
    proc1 = ""
    qualities = proc2
    proc2 = ""    
    Q = 0

    for i in range(0, len(qualities)):
        #print(str(qualities[i])+"\t"+str(ord(qualities[i])-33))
        if (ord(qualities[i])-33) > Q:
            proc1 = proc1 + bases[i]
            proc2 = proc2 + qualities[i]
    
    return proc1

#work directory is 
import argparse

parser = argparse.ArgumentParser() #make a parser

parser.add_argument("input", metavar = "prefix", default = None, type = str,
                        help = "prefix of input file") 
    
parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                        help = "Path to input file")

parser.add_argument("--genome_id", choices = ["hg19", "hg38"], default = "hg19",
                          help = "Genome id used for selecting snp data (default: %(default)s)") 

parser.add_argument("--gbam_chr_prefix", choices = ["chr", "none"], default = "none",
                          help = "chr prefix used in your bam (default: %(default)s)") 

parser.add_argument("--gbam", metavar = "wgs_bam", default = "g_bam_list.txt", type = str,
                        help = "Path to wgs bam folder.") 

args = parser.parse_args()

pr = args.input
fo = args.folder
gbamchr = args.gbam_chr_prefix

Q = 0
fo1 = "./alterativeSJ_mutprediction/"+fo+"/"
input_f0 = fo1 + pr + ".SJ.fil.annot.assadjunifreqT.pmut.SJinSJ.snp.rmut.txt"
input_file = fo1 + pr + "_g_position.txt"

if args.genome_id == "hg19" and gbamchr == "none":
    reference = "./reference/GRCh37.fa"
    #reference = "/Volumes/NaIIDA_2018aug/genome/hg19_0chr.fa"
elif args.genome_id == "hg19" and gbamchr == "chr":
    reference = "./reference/GRCh37_chr.fa"
elif args.genome_id == "hg38" and gbamchr == "none": #chr prefix is none.
    reference = "./reference/GRCh38.d1.vd1.fa"
    #reference = "/Volumes/NaIIDA_2018aug/genome/hg38_chr.fa"
elif args.genome_id == "hg38" and gbamchr == "chr": #chr prefix is chr.
    reference = "./reference/GRCh38_chr.fa"

tmp = fo1 + pr + "_g_mpileuped.txt"

if_path = pathlib.Path(input_f0)

m_file = fo1 + if_path.stem + ".g_mpileup_ori.txt" #sample position pileup info
m2_file = fo1 + if_path.stem + ".g_mpileup_tidy{}Q.txt".format(Q)
out_file = fo1 + if_path.stem + ".gmut.txt"

if os.path.exists(args.gbam):    
    data = []
    with open(args.gbam,"r") as fi: #<sample name of junction file>,<sample name of bam> 
        reader = csv.reader(fi)
        for row in reader:
            
            data.append(row)
    datadict = dict(data)
    bam = datadict.get(pr)

else:
    print("no g_bam_list.txt")


#make position file
qf = pd.read_csv(input_f0, sep='\t', header=0, index_col=None, dtype = 'object',usecols=[0,5,16,18])
qf.columns = ['CHR','sample', 'POS','MUT']
qf =qf.loc[:,['sample','CHR', 'POS','MUT']]
qf1 = qf.drop_duplicates()
qf2 = qf1.query('POS != "-"')
qf2.to_csv(input_file, index=False, sep='\t', header=False)

with open(input_file, 'r') as in1, open(m_file, 'w') as mout:
    for line in in1:
        if line == "":
            continue
        l = line.rstrip('\n')
        F = line.rstrip('\n').split('\t')
        
        if not F[1].startswith('chr') and gbamchr == "none":
            position = F[1] + ':' + F[2] +'-'+ F[2]
        elif not F[1].startswith('chr') and gbamchr == "chr":
            position = 'chr' + F[1] + ':' + F[2] +'-'+ F[2]
        elif F[1].startswith('chr') and gbamchr == "none":
            F1 = (F[1].replace('chr', ''))
            position = F1 + ':' + F[2] +'-'+ F[2]
        elif F[1].startswith('chr') and gbamchr == "chr":
            position = F[1] + ':' + F[2] +'-'+ F[2]
        print(position)
        #samtools mpileup -r "12:123873242-123873242" -f s3://niida-tokyo/lung_wg/hg19_0chr.fa "s3://niida-tokyo/lung_wg/A549.markdup.bam" > mpileup.txt
        mpileup_commands = ["samtools", "mpileup", "-r", position, "-f", reference, bam, "-O", "-o", tmp]
        #mpileup_commands = ["/usr/local/bin/samtools", "mpileup", "-r", "6:43469413-43469413", "-f", reference, bam_folder +"H2126.Aligned.sortedByCoord.out.bam"]
        subprocess.run(mpileup_commands) 
        
        with open(tmp, 'r') as in2:
            col = in2.read()
            if col == "":
                continue
            else:
                mout.write(F[0] + "\t" + str(col))

            
#arrange of mpileup file            
with open(m_file, 'r') as in3, open(m2_file, 'w') as m2out:
    for line in in3:  #sample, chr, pos, ,ref, depth, bases, Q, readsposition
        if line == "":
            continue
        col = line.rstrip('\n').split('\t')
        print(col)
        base = tidy_bases(col[5], col[6]) 
        
        depth = 0
        base2num = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 0}
        #base2pos = {'A': [], 'C': [], 'G': [], 'T': [], 'N': [], 'a': [], 'c': [], 'g': [], 't': [], 'n': []}
        rec1 = ""
        rec2 = ""
        rec3 = ""
        rec4 = ""
        for i in range(len(base)):
            depth = depth + 1
            if base[i] == '.':
                base2num[col[3].upper()] = base2num[col[3].upper()] + 1 #col[3]=REF
                #base2pos[F[3].upper()].append(pos_vector[i])
            elif base[i] == ',':
                base2num[col[3].lower()] = base2num[col[3].lower()] + 1 
                #base2pos[F[3].lower()].append(pos_vector[i])
            else:
                base2num[base[i]] = base2num[base[i]] + 1
                #base2pos[F5[i]].append(pos_vector[i])

            if depth == 0: continue
                
            nuc = col[3].upper() #REF
            depth_p = base2num['A'] + base2num['C'] + base2num['G'] + base2num['T']
            depth_n = base2num['a'] + base2num['c'] + base2num['g'] + base2num['t']
                
            A_rate = float(base2num['A'] + base2num['a'])/(depth_p + depth_n)
            A_reads = base2num['A'] + base2num['a']
            T_rate = float(base2num['T'] + base2num['t'])/(depth_p + depth_n)
            T_reads = base2num['T'] + base2num['t']
            G_rate = float(base2num['G'] + base2num['g'])/(depth_p + depth_n)
            G_reads = base2num['G'] + base2num['g']
            C_rate = float(base2num['C'] + base2num['c'])/(depth_p + depth_n)
            C_reads = base2num['C'] + base2num['c']
                
                #var_p = base2num[alt.upper()]
                #var_n = base2num[alt.lower()]
                #strand_ratio =  str(base2num[nuc]) + ":" + str(base2num[nuc.lower()])
                #rec = '\t'.join(F) + '\t' + str(depth_p + depth_n) + '\t' + \
                #str(var_p + var_n) + '\t' + str(round(alt_rate, 4)) + '\t' + str(strand_ratio)+'\n'
                
            rec1 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tA\t" + base + "\t" + str(A_reads) + "\t" + str(A_rate) + "\n"
            rec2 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tT\t" + base + "\t" + str(T_reads) + "\t" + str(T_rate) + "\n"
            rec3 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tG\t" + base + "\t" + str(G_reads) + "\t" + str(G_rate) + "\n"
            rec4 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tC\t" + base + "\t" + str(C_reads) + "\t" + str(C_rate) + "\n"
                
        m2out.write(rec1)
        m2out.write(rec2)
        m2out.write(rec3)
        m2out.write(rec4)

            
def f(row):
    if (str(row['g_bases']) != '-') & (int(row['g_alt_reads'])>1) & (float(row['g_alt_ratio'])>0):
        val = "True"
    elif str(row['g_bases']) == '-':
        val = "na"
    else:
        val = "False"
    return val

size = os.path.getsize(m2_file)

if size != 0:
    
    df1 = pd.read_csv(m2_file, sep='\t', header=None, index_col=None, dtype = 'object')
    df1.columns = ['sample','CHR', 'POS', 'MUT', 'g_bases', 'g_alt_reads', 'g_alt_ratio']
    df2 = pd.read_csv(input_f0, sep='\t', header=0, index_col=None, dtype = 'object')
    
    if not df1.iloc[0,1].startswith('chr') and not df2.iloc[0,0].startswith('chr') : #m_pileup=none, none
        print("chr_prefix:bam=none, jounction=none")
    elif not df1.iloc[0,1].startswith('chr') and df2.iloc[0,0].startswith('chr'): #none, chr
        print("chr_prefix:bam=none, jounction=chr")
        df1['pre'] = "chr"
        df1['CHR'] = df1['pre'] + df1['CHR']
        df1 = df1.drop('pre', axis=1)
        
    elif df1.iloc[0,1].startswith('chr') and not df2.iloc[0,0].startswith('chr'): #chr,none
        print("chr_prefix:bam=chr, jounction=none")
        df1['CHR'] = df1['CHR'].str.split('chr', expand=True)[1]
        
    elif df1.iloc[0,1].startswith('chr') and df2.iloc[0,0].startswith('chr'):
        print("chr_prefix:bam=chr, jounction=chr")
    
    res = pd.merge(df2, df1, on=['sample','CHR', 'POS', 'MUT'],how='left').drop_duplicates()
    res=res.fillna({'g_bases': '-', 'g_alt_reads': 0, 'g_alt_ratio': 0})
    
    res['g_mut'] = res.apply(f, axis=1)
    
    res.to_csv(out_file, index=False, sep='\t', header=True)

else:
    
    df2 = pd.read_csv(input_f0, sep='\t', header=0, index_col=None, dtype = 'object')
    df2['rna_bases'] = '-'
    df2['rna_alt_reads'] = '0'
    df2['rna_alt_ratio'] = '0'
    res = df2.drop_duplicates()
    
    res['rna_mut'] = 'na'
    
    res.to_csv(out_file, index=False, sep='\t', header=True)


os.remove(input_file)
os.remove(tmp)


