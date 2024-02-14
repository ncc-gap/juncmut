#! /usr/bin/env python3

"""
The program check bam or cram from the filename.
"""

import os
import csv
import pysam
import subprocess

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
def tidy_reads(seq, qualities, read_ids, mut_mut):
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
        if proc[i].upper() == mut_mut:
            if (ord(qualities[i])-33) > Q:
                pos_read_id_list.append(read_id_list[i])
                
    return pos_read_id_list

def juncmut_supportread_count(input_file, output_file, input_bam_file, reference):

    if input_bam_file.endswith('.bam'):
        bam_file = pysam.AlignmentFile(input_bam_file, 'rb')
    elif input_bam_file.endswith('.cram'):
        bam_file = pysam.AlignmentFile(input_bam_file, 'rc')

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Support_read_rmdup", "RNA_mut"
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:
            # Is a position of mutation in Exon or Intron
            if csvobj["Realign_result"] != "True":
                continue

            # mpileup
            (mut_chr, mut_pos, mut_ref, mut_mut) = csvobj["Mut_key"].split(',')
            mpileup_commands = ["samtools", "mpileup", "-r", "%s:%s-%s" % (mut_chr, mut_pos, mut_pos), "-f", reference, input_bam_file, "--output-QNAME", "-o", output_file + ".tmp1.txt"]
            subprocess.run(mpileup_commands)
            
            # extract read id with mutations.
            pos_read_list = []
            with open(output_file + ".tmp1.txt", 'r') as hin_tmp:
                for line in hin_tmp: 
                    col = line.rstrip('\n').split('\t')
                    bases = col[4]
                    qualities = col[5]
                    read_ids = col[6]
                    
                    reads_with_mut_list = tidy_reads(bases, qualities, read_ids, mut_mut)
                
                for read in bam_file.fetch(region = "%s:%s-%s" % (mut_chr, mut_pos, mut_pos)):
                    if not check_read(read):
                        continue
                    if read.qname in reads_with_mut_list:
                        pos_read = "%d_%d_%d" % (read.reference_start, read.reference_end, read.next_reference_start)
                        pos_read_list.append(pos_read)

            support_read_rmdup = len(set(pos_read_list))
            if support_read_rmdup >= 2:
                csvobj["Support_read_rmdup"] = support_read_rmdup
                csvobj["RNA_mut"] = "True"
                csvwriter.writerow(csvobj)

            os.remove(output_file + ".tmp1.txt")

    bam_file.close()

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    input_bam_file = sys.argv[3]
    reference = sys.argv[4]

    juncmut_supportread_count(input_file, output_file, input_bam_file, reference)
    
"""
tidy_bases(bases, qualities)
bases = "<<>><>><><><<<<>><<>>>>>>G>>><>>>><<>>>><>><<>>>>>>>>>CCC"
qualities =	"FFmcFllHDJmJJJJHsIJJiFJI7JmJGJJJk>JJJsJJIFJCDH7FFFDFDFJFF"
"""
