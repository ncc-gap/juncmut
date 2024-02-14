#! /usr/bin/env python

# [TODO] switch this function
"""
Q = 0
def tidy_bases(bases, qualities):
    output_bases = ""
    output_qualities = ""

    while bases != '':
        if bases[0] in ['>', '<', '*']: 
            bases = bases[1:]
            qualities = qualities[1:]
        elif bases[0] in '^':
            bases = bases[2:]
        elif bases[0] in '$':
            bases = bases[1:]
        elif bases[0] in ['.', ',', 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
            if (ord(qualities[0])-33) > Q:
                output_bases += bases[0]
                output_qualities += qualities[0]

            bases = bases[1:]
            qualities = qualities[1:]
            if len(bases) > 0 and bases[0] in ['+', '-']:
                match = re.search(r'^[\+\-](\d+)', bases)
                indel_size = int(match.group(1))
                bases = bases[len(str(indel_size)) + indel_size + 1:]

    if len(output_bases) != len(output_qualities):
        print("Error???")
        sys.exit(1)

    return output_bases
"""
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

def juncmut_rnamut(input_file, output_file, rna_bam, reference):
    import subprocess
    import csv
    import os

    # separate records for each variant and create position list. If no candidates, create zero size file.
    input_header = []
    data_raw = {}
    with open(input_file) as hin, open(output_file + ".tmp1.pos.bed", 'w') as hout_tmp1:
        csvreader = csv.DictReader(hin, delimiter='\t')
        input_header = csvreader.fieldnames
        for csvobj in csvreader:
            if csvobj["Possive_alt_key"] == '-':
                continue
            for alt_key in csvobj["Possive_alt_key"].split(','):
                (alt_key_pos, alt_key_ref, alt_key_alt) = alt_key.split(":")
                for var in alt_key_alt:
                    mut_key_raw = "%s,%s,%s,%s,%s,%s" % (csvobj["Chr"], csvobj["Start"], csvobj["End"], alt_key_pos, alt_key_ref, var)
                    data_raw[mut_key_raw] = csvobj
                    print('\t'.join([csvobj["Chr"], str(int(alt_key_pos) - 1), alt_key_pos, var]), file = hout_tmp1)
    # mpileup
    if data_raw != {}:
        mpileup_commands = ["samtools", "mpileup", "-l", output_file + ".tmp1.pos.bed", "-f", reference, rna_bam,"-O", "-o", output_file + ".tmp2"]
        subprocess.run(mpileup_commands)
    else:
        open(output_file + ".tmp2", "w").close()

    # arrange of mpileup file
    data_mpileup = {}
    with open(output_file + ".tmp2") as hin_tmp2:
        for line in hin_tmp2:  
            (m_chr, m_pos, m_ref, _, m_bases, m_quality, _) = line.rstrip('\n').split('\t')
            base = tidy_bases(m_bases, m_quality)

            base2num = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 0}
            A_freq = 0
            T_freq = 0
            G_freq = 0
            C_freq = 0
            A_reads = 0
            T_reads = 0
            G_reads = 0
            C_reads = 0
            for i in range(len(base)):
                if base[i] == '.':
                    base2num[m_ref.upper()] = base2num[m_ref.upper()] + 1 
                elif base[i] == ',':
                    base2num[m_ref.lower()] = base2num[m_ref.lower()] + 1
                else:
                    base2num[base[i]] = base2num[base[i]] + 1

                depth = base2num['A'] + base2num['C'] + base2num['G'] + base2num['T'] + base2num['a'] + base2num['c'] + base2num['g'] + base2num['t']

                A_freq = float(base2num['A'] + base2num['a'])/depth
                A_reads = base2num['A'] + base2num['a']
                T_freq = float(base2num['T'] + base2num['t'])/depth
                T_reads = base2num['T'] + base2num['t']
                G_freq = float(base2num['G'] + base2num['g'])/depth
                G_reads = base2num['G'] + base2num['g']
                C_freq = float(base2num['C'] + base2num['c'])/depth
                C_reads = base2num['C'] + base2num['c']

            data_mpileup["%s,%s,%s,A" % (m_chr, m_pos, m_ref)] = {"base": base, "reads": A_reads, "freq": A_freq}
            data_mpileup["%s,%s,%s,T" % (m_chr, m_pos, m_ref)] = {"base": base, "reads": T_reads, "freq": T_freq}
            data_mpileup["%s,%s,%s,G" % (m_chr, m_pos, m_ref)] = {"base": base, "reads": G_reads, "freq": G_freq}
            data_mpileup["%s,%s,%s,C" % (m_chr, m_pos, m_ref)] = {"base": base, "reads": C_reads, "freq": C_freq}

    with open(output_file, 'w') as hout:
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=input_header + ["Mut_pos", "Mut_ref", "Mut_alt", "Mut_count", "Mut_depth", "Mut_freq"])
        csvwriter.writeheader()
        for mut_key_raw in data_raw:
            output_obj = {}
            for col_name in data_raw[mut_key_raw]:
                output_obj[col_name] = data_raw[mut_key_raw][col_name]

            chrom, _, _, alt_key_pos, alt_key_ref, var = mut_key_raw.split(",")
            mut_key = "%s,%s,%s,%s" % (chrom, alt_key_pos, alt_key_ref, var)
            output_obj["Mut_pos"] = alt_key_pos
            output_obj["Mut_ref"] = alt_key_ref
            output_obj["Mut_alt"] = var
            #output_obj["Mut_bases"] = '-'
            # [TODO] fix to 0
            output_obj["Mut_depth"] = 1
            output_obj["Mut_count"] = 0
            output_obj["Mut_freq"] = 0
            if mut_key in data_mpileup:
                #output_obj["Mut_bases"] = data_mpileup[mut_key]["base"]
                # [TODO] remove this code
                base = data_mpileup[mut_key]["base"]
                if base == "":
                    base = "-"
                output_obj["Mut_depth"] = len(base)
                output_obj["Mut_count"] = data_mpileup[mut_key]["reads"]
                output_obj["Mut_freq"] = data_mpileup[mut_key]["freq"]
            csvwriter.writerow(output_obj)

    os.remove(output_file + ".tmp1.pos.bed")
    os.remove(output_file + ".tmp2")

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bam_file = sys.argv[3]
    reference = sys.argv[4]

    juncmut_rnamut(input_file, output_file, bam_file, reference)

