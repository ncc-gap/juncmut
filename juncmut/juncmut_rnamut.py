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
# -------------------------------------------
# juncmut_rnamut's function
# -------------------------------------------
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

# -------------------------------------------
# juncmut_supportread_count's function
# -------------------------------------------
def check_read(read):
    check_flag = True 

    # get the flag information
    flags = format(int(read.flag), "#014b")[:1:-1]

    # skip unmapped read  
    if flags[2] == "1" or flags[3] == "1":
        check_flag = False

    # skip supplementary alignmentx
    elif flags[8] == "1" or flags[11] == "1":
        check_flag = False

    # skip duplicated reads
    elif flags[10] == "1":
        check_flag = False

    return check_flag

# [TODO] use common module
def tidy_reads(seq, qualities, read_ids, mut):
    import re
    read_id_list = read_ids.split(',')
    proc = ""
    Q = 15
    pos_read_id_list = []
    seq_length = len(seq)
    baseIndex = 0
    # modify seq to 1 base presented by a char.
    while baseIndex < seq_length:
        #A ’>’ or ’<’ for a reference skip.
        # The deleted bases will be presented as ‘*’ in the following lines. 
        if seq[baseIndex] == '>' or seq[baseIndex] == '<' or seq[baseIndex] == '*' :
           proc = proc + seq[baseIndex]
           baseIndex += 1  
        #A '^' the end or start of read, following the quality and the base.
        elif seq[baseIndex] == '^':
            proc = proc + seq[baseIndex+2]
            baseIndex += 3
        #A '$' is the last position of read. 
        elif seq[baseIndex] == '$':
            baseIndex += 1
        #\+[0-9]+[bases] or -[0-9]+[bases] means the deletion and the insertion. For example, +2AG means insertion of AG in the forward strand
        elif seq[baseIndex] == '+' or seq[baseIndex] == '-':
            indel_length = re.search(r'\d+', seq[baseIndex:]).group()
            baseIndex += len(str(indel_length))+int(indel_length)+1 
        else:
            proc = proc + seq[baseIndex]
            baseIndex += 1
            
    # quality and base check. extract id.
    for i in range(0, len(proc),1):
        if proc[i].upper() == mut:
            if (ord(qualities[i])-33) > Q:
                pos_read_id_list.append(read_id_list[i])
                
    return pos_read_id_list

def juncmut_rnamut(input_file, output_file, input_bam_file, reference, mut_num_thres, mut_freq_thres):
    import subprocess
    import csv
    import os
    import pysam

    # separate records for each variant and create position list. If no candidates, create zero size file.
    input_header = []
    mut_key_list = []
    with open(input_file) as hin, open(output_file + ".tmp1.pos.bed", 'w') as hout_tmp1:
        csvreader = csv.DictReader(hin, delimiter='\t')
        input_header = csvreader.fieldnames
        for csvobj in csvreader:
            if csvobj["Possive_alt_key"] == '-':
                continue
            (sj_key_chr, sj_key_start, sj_key_end) = csvobj["SJ_key"].split(",")
            for alt_key in csvobj["Possive_alt_key"].split(','):
                (alt_key_pos, alt_key_ref, alt_key_alt) = alt_key.split(":")
                for var in alt_key_alt:
                    hout_tmp1.write("%s\t%d\t%s\t%s\n" % (sj_key_chr, int(alt_key_pos) - 1, alt_key_pos, var))
                    mut_key_list.append("%s,%s,%s,%s" % (sj_key_chr, alt_key_pos, alt_key_ref, var))

    # mpileup
    if len(mut_key_list) > 0:
        subprocess.run(["samtools", "mpileup", "-l", output_file + ".tmp1.pos.bed", "-f", reference, input_bam_file, "--output-QNAME", "-o", output_file + ".tmp2"])
    else:
        open(output_file + ".tmp2", "w").close()

    if input_bam_file.endswith('.bam'):
        bam_file = pysam.AlignmentFile(input_bam_file, 'rb')
    elif input_bam_file.endswith('.cram'):
        bam_file = pysam.AlignmentFile(input_bam_file, 'rc')

    # arrange of mpileup file
    data_mpileup = {}
    with open(output_file + ".tmp2") as hin_tmp2:
        for line in hin_tmp2:  
            (mpileup_chr, mpileup_pos, mpileup_ref, _, mpileup_bases, mpileup_quality, mpileup_read_ids) = line.rstrip('\n').split('\t')
            base = tidy_bases(mpileup_bases, mpileup_quality)

            base2num = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
            for i in range(len(base)):
                if base[i] in ['.', ',']:
                    base2num[mpileup_ref.upper()] += 1
                else:
                    base2num[base[i].upper()] += 1

            for mut_alt in ['A', 'C', 'G', 'T']:
                mut_key = "%s,%s,%s,%s" % (mpileup_chr, mpileup_pos, mpileup_ref, mut_alt)
                if not mut_key in mut_key_list:
                    continue

                # ↓juncmut_supportread_count
                reads_with_mut_list = tidy_reads(mpileup_bases, mpileup_quality, mpileup_read_ids, mut_alt)
                pos_read_list = []
                for read in bam_file.fetch(region = "%s:%s-%s" % (mpileup_chr, mpileup_pos, mpileup_pos)):
                    if not check_read(read):
                        continue
                    if read.qname in reads_with_mut_list:
                        pos_read = "%d_%d_%d" % (read.reference_start, read.reference_end, read.next_reference_start)
                        pos_read_list.append(pos_read)
                data_mpileup[mut_key] = {"base": base, "reads": base2num[mut_alt], "support_read_rmdup": len(set(pos_read_list))}
    bam_file.close()

    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + ["Mut_key", "Pileup_mut_count", "Pileup_depth", "Support_read_rmdup"])
        csvwriter.writeheader()

        for csvobj in csvreader:
            if csvobj["Possive_alt_key"] == '-':
                continue
            (sj_key_chr, sj_key_start, sj_key_end) = csvobj["SJ_key"].split(",")
            for alt_key in csvobj["Possive_alt_key"].split(','):
                (alt_key_pos, alt_key_ref, alt_key_alt) = alt_key.split(":")
                for var in alt_key_alt:
                    mut_key = "%s,%s,%s,%s" % (sj_key_chr, alt_key_pos, alt_key_ref, var)
                    if not mut_key in data_mpileup:
                        continue

                    # [TODO] remove this code, fix to 0
                    base = data_mpileup[mut_key]["base"]
                    if base == "":
                        base = "-"

                    pileup_depth = len(base)
                    pileup_mut_count = data_mpileup[mut_key]["reads"]
                    support_read_rmdup = data_mpileup[mut_key]["support_read_rmdup"]

                    #if RNA_mutation reads>=2 and Freq>=0.05, do realign.
                    if float(pileup_mut_count/pileup_depth) < mut_freq_thres or pileup_mut_count < mut_num_thres or support_read_rmdup < 2:
                        continue

                    csvobj["Mut_key"] = mut_key
                    csvobj["Pileup_depth"] = pileup_depth
                    csvobj["Pileup_mut_count"] = pileup_mut_count
                    csvobj["Support_read_rmdup"] = support_read_rmdup

                    csvwriter.writerow(csvobj)

    os.remove(output_file + ".tmp1.pos.bed")
    os.remove(output_file + ".tmp2")

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bam_file = sys.argv[3]
    reference = sys.argv[4]
    mut_num_thres = 2
    mut_freq_thres = 0.05 
    juncmut_rnamut(input_file, output_file, bam_file, reference, mut_num_thres, mut_freq_thres)

