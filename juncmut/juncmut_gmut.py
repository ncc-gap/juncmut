#! /usr/bin/env python
"""
Naoko Iida

python juncmut_annotgmut.py input_file, output_file, bam, reference, is_grc
if is_grc=="True", remove "chr" in input.
if is_grc=="False", add "chr" in output.
"""

def juncmut_gmut(input_file, output_file, dna_bam, reference, is_grc, mut_num_thres, mut_freq_thres):  

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
            
            if (ord(qualities[i])-33) > Q:
                proc1 = proc1 + bases[i]
                proc2 = proc2 + qualities[i]
        
        return proc1

    ##mpileup

    # separate records for each variant and create position list
    with open(input_file, 'r') as hin, open(output_file + ".tmp1", 'w') as hout1, open(output_file + ".tmp2", 'w') as hout2:
        header = hin.readline()
        header_list = header.rstrip('\n').split('\t')
        ncol = len(header_list)
        print(ncol)
        header_list[0] = "Mut_key"
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mut_key = F[0]
            P = mut_key.split(',')
            mut_pos = P[1]
           
            if is_grc == 'True':
                chr = P[0].replace('chr', '')
            else:
                chr_t = P[0].replace('chr', '')
                chr = 'chr'+chr_t
                
            new_mut_key = chr+","+str(mut_pos)+","+P[2]+","+P[3]
            
            hout1.write(new_mut_key+'\t'+'\t'.join(F[1:int(ncol)])+'\n')
            
            position = str(chr) + ':' + str(int(mut_pos)-1) + '-' + str(mut_pos)
            #import pdb; pdb.set_trace()
            mpileup_commands = ["samtools", "mpileup", "-r", position, "-f", reference, dna_bam, "-O", "-o", output_file + ".tmp22"]
            subprocess.run(mpileup_commands) 
        
            with open(output_file + ".tmp22", 'r') as m2tmp:
                col = m2tmp.read()
                if col == "":
                    continue
                else:
                    hout2.write(str(col))


    #arrange of mpileup file
    with open(output_file + ".tmp2", 'r') as in3, open(output_file + ".tmp3", 'w') as m2out:
        for line in in3:
            
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

                #rec1 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tA\t" + base + "\t" + str(A_reads) + "\t" + str(A_rate) + "\n"
                #rec2 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tT\t" + base + "\t" + str(T_reads) + "\t" + str(T_rate) + "\n"
                #rec3 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tG\t" + base + "\t" + str(G_reads) + "\t" + str(G_rate) + "\n"
                #rec4 = col[0] +"\t"+ col[1] +"\t"+ col[2] + "\tC\t" + base + "\t" + str(C_reads) + "\t" + str(C_rate) + "\n"
                
                rec1 = col[0] +","+ col[1] +","+ col[2] + ",A\t" + str(len(base)) + "\t" + str(A_reads) + "\t" + str(A_rate) + "\n"
                rec2 = col[0] +","+ col[1] +","+ col[2] + ",T\t" + str(len(base)) + "\t" + str(T_reads) + "\t" + str(T_rate) + "\n"
                rec3 = col[0] +","+ col[1] +","+ col[2] + ",G\t" + str(len(base)) + "\t" + str(G_reads) + "\t" + str(G_rate) + "\n"
                rec4 = col[0] +","+ col[1] +","+ col[2] + ",C\t" + str(len(base)) + "\t" + str(C_reads) + "\t" + str(C_rate) + "\n"

            m2out.write(rec1)
            m2out.write(rec2)
            m2out.write(rec3)
            m2out.write(rec4)


    def f(row):
        if (int(row['g_depth']) > 0) & (int(row['g_alt_reads']) >= mut_num_thres) & (float(row['g_alt_ratio']) >= mut_freq_thres):
            val = "True"
        elif int(row['g_depth']) == 0:
            val = "False"
        else:
            val = "False"
        return val

    fsize = Path(output_file + ".tmp3").stat().st_size

    if fsize != 0:
        df1 = pd.read_csv(output_file+".tmp3", sep='\t', header=None, index_col=None, dtype = 'object')
        df1.columns = ['Mut_key_', 'g_depth', 'g_alt_reads', 'g_alt_ratio']
        
        df2 = pd.read_csv(output_file+".tmp1", sep='\t', header=None, index_col=None, dtype = 'object')
        
        df2.columns = header_list
        
        res = pd.merge(df2, df1, on=['Mut_key_'],how='left').drop_duplicates()
        res=res.fillna({'g_depth': 0, 'g_alt_reads': 0, 'g_alt_ratio': 0})

        res['g_mut'] = res.apply(f, axis=1)

        res.to_csv(output_file, index=False, sep='\t')

    else:

        df2 = pd.read_csv(output_file + ".tmp1", sep='\t', header=None, index_col=None, dtype = 'object')
        df2.columns = header_list
        df2['g_depth'] = '-'
        df2['g_alt_reads'] = '0'
        df2['g_alt_ratio'] = '0'
        res = df2.drop_duplicates()

        res['g_mut'] = 'False'
        res.to_csv(output_file, index=False, sep='\t', header=True)

    Path(output_file + ".tmp1").unlink()
    Path(output_file + ".tmp2").unlink()
    Path(output_file + ".tmp3").unlink()
    Path(output_file + ".tmp22").unlink()

    
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
    
    parser.add_argument("--is_grc", metavar = "is_grc", default = "True", type = str,
                            help = "True means no chr prefix in bam. False means chr prefix in bam")
    
    parser.add_argument("--mut_num_thres", type = int, default = 2,
                        help = "A mutation with mutation alleles >= mut_num_thres is a true candidate (default: %(default)s)")
    
    parser.add_argument("--mut_freq_thres", type = float, default = 0.05,
                        help = "A mutation with frequency >= mut_freq_thres is a true candidate (default: %(default)s)")

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    dna_bam = args.dna_bam
    reference = args.reference
    is_grc = args.is_grc
    mut_num_thres = args.mut_num_thres 
    mut_freq_thres = args.mut_freq_thres

    juncmut_gmut(input_file, output_file, dna_bam, reference, is_grc, mut_num_thres, mut_freq_thres)
