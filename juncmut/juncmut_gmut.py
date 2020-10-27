#! /usr/bin/env python
"""
Naoko Iida

python juncmut_annotgmut.py input_file, output_file, bam, reference, is_grc
if is_grc=="T", remove "chr" in input.
if is_grc=="F", add "chr" in output.
"""

def juncmut_gmut(input_file, output_file, dna_bam, reference, is_grc):
    
    import subprocess
    import pandas as pd
    from pathlib import Path

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
        Q = 0
    
        for i in range(0, len(qualities)):
            #print(str(qualities[i])+"\t"+str(ord(qualities[i])-33))
            if (ord(qualities[i])-33) > Q:
                proc1 = proc1 + bases[i]
                proc2 = proc2 + qualities[i]
        
        return proc1

    ##mpileup

    # separate records for each variant and create position list
    with open(input_file, 'r') as hin, open(output_file + ".tmp1", 'w') as hout1, open(output_file + ".tmp1.pos.bed", 'w') as hout2:
        next(hin)
        for line in hin:            
            F = line.rstrip('\n').split('\t')
            mut_pos = F[14]
            var = F[16]           
            if is_grc == 'T':
                chr = F[0].replace('chr', '')
            else:
                chr_t = F[0].replace('chr', '')
                chr = 'chr'+chr_t
            print(str(chr)+'\t'+'\t'.join(F[1:3]), file = hout1)
            print('\t'.join([chr, str(int(mut_pos) - 1), str(mut_pos), var]), file = hout2)
            
    mpileup_commands = ["samtools", "mpileup", "-l", output_file + ".tmp1.pos.bed", "-f", reference, dna_bam, "-O", "-o", output_file + ".tmp2"]      
    subprocess.run(mpileup_commands) 

    #arrange of mpileup file
    with open(output_file + ".tmp2", 'r') as in3, open(output_file + ".tmp3", 'w') as m2out:
        for line in in3:  #sample, chr, pos, ,ref, depth, bases, Q, readsposition
            
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
                    base2num[col[2].upper()] = base2num[col[2].upper()] + 1 #col[3]=REF
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


    def f(row):
        if (str(row['g_bases']) != '-') & (int(row['g_alt_reads'])>1) & (float(row['g_alt_ratio'])>0.05):
            val = "T"
        elif str(row['g_bases']) == '-':
            val = "na"
        else:
            val = "F"
        return val


    # size = os.path.getsize(output_file + ".tmp3")
    fsize = Path(output_file + ".tmp3").stat().st_size

    if fsize != 0:
        df1 = pd.read_csv(output_file+".tmp3", sep='\t', header=None, index_col=None, dtype = 'object')
        df1.columns = ['Chr', 'Mut_Pos', 'Mut_Ref', 'Mut_Alt', 'g_bases', 'g_alt_reads', 'g_alt_ratio']
        
        df2 = pd.read_csv(output_file+".tmp1", sep='\t', header=None, index_col=None, dtype = 'object')
        
        df2.columns = ['Chr','SJ_Start','SJ_End','SJ_Type','SJ_Strand','SJ_Read_Count','SJ_Depth','SJ_Freq','Ref_Motif','Possivle_Alt_Motif',
'Possible_Alt_key','Is_Canonical','Is_in_exon','SJ_Overlap_Count','Mut_Pos','Mut_Ref','Mut_Alt','Mut_Count','Mut_Depth','Mut_Freq',
'Realign_No_SJ_Neg','Realign_No_SJ_Pos','Realign_Target_SJ_Neg','Reaglin_Target_SJ_Pos','Realign_Normal_SJ_Neg','Realign_Normal_SJ_Pos','RNA_Mut','gnomAD','gnomAD_AF']
        
        res = pd.merge(df2, df1, on=['Chr', 'Mut_Pos', 'Mut_Ref', 'Mut_Alt'],how='left').drop_duplicates()
        res=res.fillna({'g_bases': '-', 'g_alt_reads': 0, 'g_alt_ratio': 0})

        res['g_mut'] = res.apply(f, axis=1)

        res.to_csv(output_file, index=False, sep='\t')

    else:

        df2 = pd.read_csv(output_file + ".tmp1", sep='\t', header=None, index_col=None, dtype = 'object')
        df2.columns = ['Chr','SJ_Start','SJ_End','SJ_Type','SJ_Strand','SJ_Read_Count','SJ_Depth','SJ_Freq','Ref_Motif','Possivle_Alt_Motif',
'Possible_Alt_key','Is_Canonical','Is_in_exon','SJ_Overlap_Count','Mut_Pos','Mut_Ref','Mut_Alt','Mut_Count','Mut_Depth','Mut_Freq',
'Realign_No_SJ_Neg','Realign_No_SJ_Pos','Realign_Target_SJ_Neg','Reaglin_Target_SJ_Pos','Realign_Normal_SJ_Neg','Realign_Normal_SJ_Pos','RNA_Mut','gnomAD','gnomAD_AF']
        df2['g_bases'] = '-'
        df2['g_alt_reads'] = '0'
        df2['g_alt_ratio'] = '0'
        res = df2.drop_duplicates()

        res['g_mut'] = 'na'
        res.to_csv(output_file, index=False, sep='\t', header=True)

    Path(output_file + ".tmp1").unlink()
    Path(output_file + ".tmp2").unlink()
    Path(output_file + ".tmp3").unlink()
    Path(output_file + ".tmp1.pos.bed").unlink()

    
if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser

    parser.add_argument("--input_file", metavar = "input_file", default = None, type = str,
                            help = "input file")

    parser.add_argument("--output_file", metavar = "output_file", default = None, type = str,
                            help = "output file")

    parser.add_argument("--dna_bam", metavar = "dna_bam", default = None, type = str,
                            help = "genomic bam")

    parser.add_argument("--reference", metavar = "reference", default = None, type = str,
                            help = "/path/to/reference")
    
    parser.add_argument("--is_grc", metavar = "is_grc", default = "T", type = str,
                            help = "T means no chr prefix in bam. F means chr prefix in bam")

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    dna_bam = args.dna_bam
    reference = args.reference
    is_grc = args.is_grc

    juncmut_gmut(input_file, output_file, dna_bam, reference, is_grc)