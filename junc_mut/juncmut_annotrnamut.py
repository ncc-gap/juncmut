#! /usr/bin/env python
"""
Naoko Iida

#chromosome prefix of junction and bam should be the same.
#sample name list is "rna_bam_list.txt".
"""

def juncmut_annotrnamut(pr, folder, genome_id, rbamchr, rbam):

    import subprocess
    import pandas as pd
    import os
    import pathlib
    import csv
    
    def tidy_bases(bases, qualities):  #remove indel and edeage
        
        import re
    
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
        Q = 15
    
        for i in range(0, len(qualities)):
            #print(str(qualities[i])+"\t"+str(ord(qualities[i])-33))
            if (ord(qualities[i])-33) > Q:
                proc1 = proc1 + bases[i]
                proc2 = proc2 + qualities[i]
        
        return proc1
    
    cdir = './data/%s/' %(folder)
    os.chdir(cdir)
    ##mpileup
    Q = 15
    fo1 = "./alterativeSJ_mutprediction/" %(folder)
    fo2 = "./alterativeSJ_mutprediction/" %(folder)
    #input_f0 = fo1 + pr + ".SJ.fil.annot.assadjunifreqT.pmut.SJinSJ.snp.txt" #original
    input_f0 = fo1 + pr + ".SJ.fil.annot.assadjunifreqT.pmut.SJinSJ.snp.txt"
    input_file = fo2 + pr + "_r_position.txt"
    
    if genome_id == "hg19" and rbamchr == "none":
        reference = "./reference/GRCh37.fa"
        #reference = "/Volumes/NaIIDA_2018aug/genome/hg19_0chr.fa"
    elif genome_id == "hg19" and rbamchr == "chr":
        reference = "./reference/GRCh37_chr.fa"
    elif genome_id == "hg38" and rbamchr == "none": #chr prefix is none.
        reference = "./reference/GRCh38.d1.vd1.fa"
        #reference = "/Volumes/NaIIDA_2018aug/genome/hg38_chr.fa"
    elif genome_id == "hg38" and rbamchr == "chr": #chr prefix is chr.
        reference = "./reference/GRCh38_chr.fa"
    
    
    tmp = fo2 + pr + "_r_mpileuped.txt"
    
    if_path = pathlib.Path(input_f0)
    
    m_file = fo2 + if_path.stem + ".rna_mpileup_ori.txt" #sample position pileup info
    m2_file = fo2 + if_path.stem + ".rna_mpileup_tidy{}Q.txt".format(Q)
    out_file = fo2 + if_path.stem + ".rmut.txt"
    
    if os.path.exists(rbam):    
        data = []
        with open(rbam,"r") as fi: #<sample name of junction file>,<sample name of bam> 
            reader = csv.reader(fi)
            for row in reader:
                
                data.append(row)
        datadict = dict(data)
        bam = datadict.get(pr)
    
    else:
        print("no rna_bam_list.txt")
    
    #make position file
    qf = pd.read_csv(input_f0, sep='\t', header=None, index_col=None, dtype = 'object',usecols=[0,5,16,18])
    qf.columns = ['CHR','sample', 'POS','MUT']
    qf =qf.loc[:,['sample','CHR', 'POS','MUT']]
    qf1 = qf.drop_duplicates()
    qf2 = qf1.query('POS != "-"')
    qf2.to_csv(input_file, index=False, sep='\t', header=False)
    
    
    with open(input_file, 'r') as in1, open(m_file, 'w') as mout:
        for line in in1:
            #l = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
    
            position = F[1] + ':' + F[2] +'-'+ F[2]
            print(position)
            #samtools mpileup -r "12:123873242-123873242" -f s3://niida-tokyo/lung_wg/hg19_0chr.fa "s3://niida-tokyo/lung_wg/A549.markdup.bam" > mpileup.txt
            mpileup_commands = ["samtools", "mpileup", "-r", position, "-f", reference, bam, "-O", "-o", tmp]
            #mpileup_commands = ["/usr/local/bin/samtools", "mpileup", "-r", "5:148630908-148630908", "-f", reference, bam_folder +"H1648.Aligned.sortedByCoord.out.bam", "-O", "-o", tmp]
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
        if (str(row['rna_bases']) != '-') & (int(row['rna_alt_reads'])>1) & (float(row['rna_alt_ratio'])>0.05):
            val = "True"
        elif str(row['rna_bases']) == '-':
            val = "na"
        else:
            val = "False"
        return val
    
    size = os.path.getsize(m2_file)
    
    if size != 0:
        
        df1 = pd.read_csv(m2_file, sep='\t', header=None, index_col=None, dtype = 'object')
        df1.columns = ['sample','CHR', 'POS', 'MUT', 'rna_bases', 'rna_alt_reads', 'rna_alt_ratio']
        #df1['CHR'] = df1['CHR'].str.split('chr', expand=True)[1]
        df2 = pd.read_csv(input_f0, sep='\t', header=None, index_col=None, dtype = 'object')
        df2.columns = ['CHR','start','end','start_ori','end_ori','sample','class','strand','reads', 'total', 'freq', 'ref_bases','mut_prediction_seq','info', 'type', 'SJinSJ','POS','REF','MUT','snp_allele','snp_freq']
        res = pd.merge(df2, df1, on=['sample','CHR', 'POS','MUT'],how='left').drop_duplicates()
        res=res.fillna({'rna_bases': '-', 'rna_alt_reads': 0, 'rna_alt_ratio': 0})
        
        res['rna_mut'] = res.apply(f, axis=1)
        
        res.to_csv(out_file, index=False, sep='\t', header=True)
    
    else:
        
        df2 = pd.read_csv(input_f0, sep='\t', header=None, index_col=None, dtype = 'object')
        df2.columns = ['CHR','start','end','start_ori','end_ori','sample','class','strand','reads', 'total', 'freq', 'ref_bases','mut_prediction_seq','info', 'type', 'SJinSJ','POS','REF','MUT','snp_allele','snp_freq']
        df2['rna_bases'] = '-'
        df2['rna_alt_reads'] = '0'
        df2['rna_alt_ratio'] = '0'
        res = df2.drop_duplicates()
        
        res['rna_mut'] = 'na'
        
        res.to_csv(out_file, index=False, sep='\t', header=True)
        
    os.remove(input_file)
    os.remove(tmp)

if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input", metavar = "prefix", default = None, type = str,
                            help = "prefix of input file") 
        
    parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                            help = "Folder to input files")
    
    parser.add_argument("--genome_id", choices = ["hg19", "hg38"], default = "hg19",
                              help = "Genome id used for selecting snp data (default: %(default)s)") 
    
    parser.add_argument("--rbam_chr_prefix", choices = ["chr", "none"], default = "none",
                              help = "chr prefix used in your bam (default: %(default)s)")
    
    parser.add_argument("--rbam", metavar = "RNAseq_bam_list", default = None, type = str,
                            help = "A file:list of Path to RNAseq bam folder.") 
        
    args = parser.parse_args()
    
    pr = args.input
    folder = args.folder
    genome_id = args.genome_id
    rbamchr = args.rbam_chr_prefix
    rbam = args.rbam
    
    juncmut_annotrnamut(pr, folder, genome_id, rbamchr, rbam)
    