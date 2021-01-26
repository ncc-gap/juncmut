#! /usr/bin/env python
"""
Naoko Iida

#chromosome prefix of junction and bam should be the same.
#sample name list is "rna_bam_list.txt".
"""

def juncmut_rnamut(input_file, output_file, rna_bam, reference):


    import subprocess
    import pandas as pd
    from pathlib import Path

    def tidy_bases(bases, qualities):  #remove indel and edeage

        import re

        proc1 = ""
        proc2 = ""

        # $ means the last position of read.  ^ is the start position of read.
        while len(bases) > 0:
            match = re.search(r'[$\^]', bases)
            
            if match is None:
                proc1 = proc1 + bases
                bases = ""
                proc2 = proc2 + qualities
                qualities = ""
            elif match.group() == "^":
                pos = match.start()
                proc1 = proc1 + bases[0:pos]
                bases = bases[(pos + 2):len(bases)]
                proc2 = proc2 + qualities[0:pos]
                qualities = qualities[(pos + 1): len(qualities)] 
            #match == "$"
            else:
                pos = match.start()
                proc1 = proc1 + bases[0:pos]
                bases = bases[(pos + 1):len(bases)]
                proc2 = proc2 + qualities[0:pos]
                qualities = qualities[(pos + 1): len(qualities)] 
                
        #remove indel. +/- +[0-9]+[atgcnATGCN*#]
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
                if match.group() == "+":
                    indel_length = re.search(r'\+\d+', bases).group().replace('+','')
                if match.group() == "-":
                    indel_length = re.search(r'\-\d+', bases).group().replace('-','')     
                skip_num = len(str(indel_length))+int(indel_length)+1
                
                proc1 = proc1 + bases[0:pos]
                bases = bases[pos+int(skip_num): len(bases)]
                proc2 = proc2 + qualities[0:pos]
                qualities = qualities[(pos + 1): len(qualities)]
        #remove skip-base position in a read.
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

##mpileup
    # separate records for each variant and create position list. If no candidates, create zero size file.
    with open(input_file, 'r') as hin, open(output_file + ".tmp1", 'w') as hout1, open(output_file + ".tmp1.pos.bed", 'w') as hout2:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[13] == '-': continue
            pmuts = F[13].split(',')
            for pmut in pmuts:
                pmut_elm = pmut.split(':')
                for var in pmut_elm[2]:
                    print('\t'.join(F + [pmut_elm[0], pmut_elm[1], var]), file = hout1)
                    print('\t'.join([F[0], str(int(pmut_elm[0]) - 1), pmut_elm[0], var]), file = hout2)
    
    query_size = Path(output_file + ".tmp1").stat().st_size
    if query_size != 0:
        
        mpileup_commands = ["samtools", "mpileup", "-l", output_file + ".tmp1.pos.bed", "-f", reference, rna_bam,"-O", "-o", output_file + ".tmp2"]
        subprocess.run(mpileup_commands)
        
    
        #arrange of mpileup file
        with open(output_file + ".tmp2", 'r') as in3, open(output_file + ".tmp3", 'w') as m2out:
            for line in in3:  
                #sample, chr, pos, ,ref, depth, bases, Q, readsposition
    
                col = line.rstrip('\n').split('\t')
    
                base = tidy_bases(col[4], col[5])
    
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
                        base2num[col[2].upper()] = base2num[col[2].upper()] + 1 
                        #col[3]=REF
                        #base2pos[F[3].upper()].append(pos_vector[i])
                    elif base[i] == ',':
                        base2num[col[2].lower()] = base2num[col[2].lower()] + 1
                        #base2pos[F[3].lower()].append(pos_vector[i])
                    else:
                        # print(base)
                        base2num[base[i]] = base2num[base[i]] + 1
                        #base2pos[F5[i]].append(pos_vector[i])
    
                    if depth == 0: continue
    
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
    
                    rec1 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tA\t" + base + "\t" + str(A_reads) + "\t" + str(A_rate) + "\n"
                    rec2 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tT\t" + base + "\t" + str(T_reads) + "\t" + str(T_rate) + "\n"
                    rec3 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tG\t" + base + "\t" + str(G_reads) + "\t" + str(G_rate) + "\n"
                    rec4 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tC\t" + base + "\t" + str(C_reads) + "\t" + str(C_rate) + "\n"
    
                m2out.write(rec1)
                m2out.write(rec2)
                m2out.write(rec3)
                m2out.write(rec4)

    else:
        
        file3 = Path(output_file + ".tmp3")
        file3.touch()
        file2 = Path(output_file + ".tmp2")
        file2.touch()
        
    # size = os.path.getsize(output_file + ".tmp3")
    fsize = Path(output_file + ".tmp3").stat().st_size

    if fsize != 0:


        df1 = pd.read_csv(output_file+".tmp3", sep='\t', header=None, index_col=None, dtype = 'object')
        df1.columns = ['CHR', 'POS', 'REF', 'MUT', 'rna_bases', 'rna_alt_reads', 'rna_alt_ratio']
        df2 = pd.read_csv(output_file+".tmp1", sep='\t', header=None, index_col=None, dtype = 'object')
        df2.columns = ['CHR','start','end','start_ori','end_ori','sample','class','strand','reads', 'total', 'freq', 'ref_bases','mut_prediction_seq','info', 'type', 'intSJ','POS','REF','MUT']
        res = pd.merge(df2, df1, on=['CHR', 'POS','REF','MUT'],how='left').drop_duplicates()
        res=res.fillna({'rna_bases': '-', 'rna_alt_reads': 0, 'rna_alt_ratio': 0})
        res.to_csv(output_file, index=False, sep='\t')

    else:
        if query_size != 0:
            df2 = pd.read_csv(output_file + ".tmp1", sep='\t', header=None, index_col=None, dtype = 'object')
            df2.columns = ['CHR','start','end','start_ori','end_ori','sample','class','strand','reads', 'total', 'freq', 'ref_bases','mut_prediction_seq','info', 'type', 'intSJ','POS','REF','MUT']
            df2['rna_bases'] = '-'
            df2['rna_alt_reads'] = '0'
            df2['rna_alt_ratio'] = '0'
            res = df2.drop_duplicates()
            res.to_csv(output_file, index=False, sep='\t')
        else:
            res = pd.DataFrame(columns=['CHR','start','end','start_ori','end_ori','sample','class','strand','reads', 'total', 'freq', 'ref_bases','mut_prediction_seq','info', 'type', 'intSJ','POS','REF','MUT','rna_bases', 'rna_alt_reads', 'rna_alt_ratio'])
            res.to_csv(output_file, index=False, sep='\t')
    Path(output_file + ".tmp1").unlink()
    Path(output_file + ".tmp2").unlink()
    Path(output_file + ".tmp3").unlink()
    Path(output_file + ".tmp1.pos.bed").unlink()

    
if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser

    parser.add_argument("-input_file", metavar = "input_file", default = None, type = str,
                            help = "input file")

    parser.add_argument("-output_file", metavar = "output_file", default = "my_sample", type = str,
                            help = "output file")

    parser.add_argument("-rna_bam", metavar = "rna_bam", default = None, type = str,
                            help = "rna bam")

    parser.add_argument("-reference", metavar = "reference", default = None, type = str,
                            help = "/path/to/reference")


    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    rna_bam = args.rna_bam
    reference = args.reference


    juncmut_rnamut(input_file, output_file, rna_bam, reference)
