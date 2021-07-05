#! /usr/bin/env python3

"""
The program check bam or cram from the filename.
"""

from pathlib import Path
import pysam
import subprocess


def juncmut_supportread_count(input_file, output_file, bam_file, reference):


    def check_read(read):

        check_flag = True 

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read  
        if flags[2] == "1" or flags[3] == "1": check_flag = False
 
        # skip supplementary alignmentx
        if flags[8] == "1" or flags[11] == "1": check_flag = False

        # skip duplicated reads
        if flags[10] == "1": check_flag = False

        return(check_flag)

   
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

    b_path = Path(bam_file)

    if b_path.suffix == '.bam':
        bamfile = pysam.AlignmentFile(bam_file, 'rb')
    if b_path.suffix == '.cram':
        bamfile = pysam.AlignmentFile(bam_file, 'rc')
 
    ## start 
    hout = open(output_file, 'w') 
    header = ["Mut_key", "SJ_key", "Sample", "SJ_Type", "SJ_Strand", "SJ_Read_Count", "SJ_Depth", "SJ_Freq",
              "Ref_Motif", "Possivle_Alt_Motif","Possible_Alt_key", "Is_GT/AG", "Is_in_exon","SJ_Overlap_Count", 
              "Chr","Mut_Pos", "Mut_Ref", "Mut_Alt", "Mut_Count", "Mut_Depth", "Mut_Freq",
              "Realign_No_SJ_Neg", "Realign_No_SJ_Pos", "Realign_Target_SJ_Neg", "Reaglin_Target_SJ_Pos",
              "Realign_Normal_SJ_Neg", "Realign_Normal_SJ_Pos","Realign_result","support_read_rmdup","RNA_Mut"]
    print('\t'.join(header), file = hout)
    # for each row.
    with open(input_file, 'r') as hin:
        next(hin)
        for line in hin:
            lie = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
            # Is a position of mutation in Exon or Intron
            if F[-1] != "True": continue
                #print(lie + "\t0\tFalse", file = hout)                
            else:
                #mpileup
                mut_elm = F[0].split(',')
                mut_chr = mut_elm[0]
                mut_pos = str(mut_elm[1])
                mut_mut = mut_elm[3]
                #samtools mpileup -r chr4:162087015-162087015 -f /Volumes/NIIDA_SSD1R/genome DRR016694.Aligned.sortedByCoord.out.bam
                mpileup_commands = ["samtools", "mpileup", "-r", mut_chr+":"+mut_pos+"-"+mut_pos, "-f", reference, bam_file, "--output-QNAME", "-o", output_file + ".tmp1.txt"]
                subprocess.run(mpileup_commands)
        
                # extract read id with mutations.
                pos_read_list = []
                with open(output_file + ".tmp1.txt", 'r') as tin:
                    for line in tin: 
                        col = line.rstrip('\n').split('\t')
                        bases = col[4]
                        qualities = col[5]
                        read_ids = col[6]
                        
                        reads_with_mut_list = tidy_reads(bases, qualities, read_ids, mut_mut)
                   
                    for read in bamfile.fetch(region = str(mut_chr) + ':' + str(mut_pos) + '-' + str(mut_pos)):
                        if not check_read(read): continue
                        else:
                            if read.qname in reads_with_mut_list:
                                pos_read = str(read.reference_start) + '_' + str(read.reference_end) + '_' + str(read.next_reference_start)   
                                pos_read_list.append(pos_read)
                support_read_rmdup = len(set(pos_read_list))
                if support_read_rmdup >= 2:
                    rna_mut = "True"
                    print(lie + "\t"+ str(support_read_rmdup) + "\t" + str(rna_mut), file = hout) 
                #else: rna_mut = "False"
                
                #print(lie + "\t"+ str(support_read_rmdup) + "\t" + str(rna_mut), file = hout) 
                
                Path(output_file + ".tmp1.txt").unlink()

    bamfile.close()
    hout.close()

    

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("-input_file", metavar = "input_file", default = None, type = str,
                            help = "input file")         
    parser.add_argument("-output_file", metavar = "output_file", default = None, type = str,
                            help = "output files") 
    parser.add_argument("-bam_file", metavar = "bam_file", default = None, type = str,
                            help = "output files") 
    parser.add_argument("-reference", metavar = "reference", default = None, type = str,
                            help = "reference")        
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    bam_file = args.bam_file
    reference = args.reference
 
    juncmut_supportread_count(input_file, output_file, bam_file, reference)
    
"""
tidy_bases(bases, qualities)
bases = "<<>><>><><><<<<>><<>>>>>>G>>><>>>><<>>>><>><<>>>>>>>>>CCC"
qualities =	"FFmcFllHDJmJJJJHsIJJiFJI7JmJGJJJk>JJJsJJIFJCDH7FFFDFDFJFF"
"""


