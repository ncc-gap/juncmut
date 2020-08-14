#! /usr/bin/env python3

"""
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/homeB/naiida/project2/juncmut_codes

The program check bam or cram from the filename.
"""

import sys, random
from pathlib import Path
import pysam
import annot_utils
from pyssw import realign_main

def juncmut_realign(input_file, output_file, bam_file, reference, genome_id, is_grc, template_size = 10, score_margin = 4):


    def generate_template_seq(template_file, ref_tb, junc_tb, mut_chr, mut_pos, mut_ref, mut_alt, 
                              junc_start, junc_end, junc_annotated, template_size, genome_id, is_grc):

        def change_base_check(tseq, tpos, tref, talt):
            # for debugging
            tseq_alt = tseq
            if tseq_alt[tpos] != tref:
                print("Inconsistent ref base", file = sys.stderr)
                sys.exit(1)
            tseq_alt = tseq[:tpos] + talt + tseq[(tpos + 1):]
            return(tseq_alt)

        hout = open(template_file, 'w') 

        ##########
        # no splicing
        # without mutation
        ir_ref_seq = ref_tb.fetch(region = mut_chr + ':' + str(mut_pos - template_size) + '-' + str(mut_pos + template_size - 1))
        print(">no_splicing_negative", file = hout)
        print(ir_ref_seq, file = hout)

        # with mutation
        ir_alt_seq = change_base_check(ir_ref_seq, template_size, mut_ref, mut_alt)
        print(">no_splicing_positive", file = hout)
        print(ir_alt_seq, file = hout)
        ##########

        ##########
        # target splicing
        # without mutation
        target_junc_ref_seq = ref_tb.fetch(region = mut_chr + ':' + str((junc_start - 1) - (template_size - 1)) + '-' + str(junc_start - 1)) + \
            ref_tb.fetch(region = mut_chr + ':' + str(junc_end + 1) + '-' + str(junc_end + 1 + (template_size - 1))) 
        print(">target_splicing_negative", file = hout)
        print(target_junc_ref_seq, file = hout)

        # with mutation
        if mut_pos >= (junc_start - 1) - (template_size - 1) and mut_pos <= junc_start - 1:
            target_junc_alt_seq = change_base_check(target_junc_ref_seq, mut_pos - ((junc_start - 1) - (template_size - 1)), mut_ref, mut_alt)
            print(">target_splicing_positive", file = hout)
            print(target_junc_alt_seq, file = hout)
        elif mut_pos >= junc_end + 1 and mut_pos <= (junc_end + 1) + (template_size - 1):
            target_junc_alt_seq = change_base_check(target_junc_ref_seq, mut_pos - ((junc_end + 1) + (template_size - 1)) - 1, mut_ref, mut_alt)
            print(">target_splicing_positive", file = hout)
            print(target_junc_alt_seq, file = hout)
        ##########

        ##########
        # normal splicing
        records = junc_tb.fetch(region = mut_chr + ':' + str(mut_pos - template_size - 1) + '-' + str(mut_pos + template_size))
        njind = 0
        for record_line in records:
            record = record_line.split('\t')
            rj_start, rj_end = int(record[1]) - 1, int(record[2])
            if (junc_annotated == junc_start and rj_start == junc_start and abs(rj_end - mut_pos) <= template_size) or \
                (junc_annotated == junc_end and rj_end == junc_end and abs(rj_start - 1 - mut_pos) <= template_size):
                # without mutation
                normal_junc_ref_seq = ref_tb.fetch(region = mut_chr + ':' + str((rj_start - 1) - (template_size - 1)) + '-' + str(rj_start - 1)) + \
                    ref_tb.fetch(region = mut_chr + ':' + str(rj_end + 1) + '-' + str((rj_end + 1) + (template_size - 1)))
                print(">normal_splicing_negative_" + str(njind), file = hout)
                print(normal_junc_ref_seq, file = hout)
            
                if mut_pos >= (rj_start - 1) - (template_size - 1) and mut_pos <= rj_start - 1:
                    normal_junc_alt_seq = change_base_check(normal_junc_ref_seq, mut_pos - ((rj_start - 1) - (template_size - 1)), mut_ref, mut_alt)
                    print(">normal_splicing_negative_" + str(njind), file = hout)
                    print(normal_junc_alt_seq)
                elif mut_pos >= rj_end + 1 and mut_pos <= (rj_end + 1) + (template_size - 1):
                    normal_junc_alt_seq = change_base_check(normal_junc_alt_seq, mut_pos - ((rj_end + 1) + (template_size - 1)) - 1, mut_ref, mut_alt)
                    print(">normal_junction_positive_" + str(njind), file = hout)
                    print(normal_junc_alt_seq, file = hout)
        ##########

        hout.close()


    def extract_read_around_boundary(bam_file, output_file, mut_chr, mut_pos, max_count = 10000):

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

        b_path = Path(bam_file)
        if b_path.suffix == '.bam':
            bamfile = pysam.AlignmentFile(bam_file, 'rb')
        if b_path.suffix == '.cram':
            bamfile = pysam.AlignmentFile(bam_file, 'rc')

        read_count = bamfile.count(region = mut_chr + ':' + str(mut_pos) + '-' + str(mut_pos), read_callback = check_read)
        read_ind_list = random.sample(range(read_count), min(max_count, read_count))
        read_inds = {}
        for i in read_ind_list:
            read_inds[i] = 1

        ind = -1
        hout = open(output_file, 'w') 
        for read in bamfile.fetch(region = mut_chr + ':' + str(mut_pos) + '-' + str(mut_pos)):

            if not check_read(read): continue

            ind = ind + 1
            if not ind in read_inds: continue

            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            read_id = read.qname + '_1' if flags[6] == '1' else read.qname + '_2'
            read_id = read_id + '_' + str(read.reference_start + 1) + '_' + str(read.reference_end) + '_' + str(read.cigarstring) 
            print('>' + read_id + '\n' + read.seq, file = hout)

        bamfile.close()
        hout.close()


    ref_tb = pysam.FastaFile(reference)
    annot_utils.junction.make_junc_info(output_file + ".gencode.junc.bed.gz", "gencode", genome_id, is_grc, False)
    junc_tb = pysam.TabixFile(output_file + ".gencode.junc.bed.gz")

    hout = open(output_file, 'w') 
    header = ["Chr", "SJ_Start", "SJ_End", "SJ_Type", "SJ_Strand", "SJ_Read_Count", "SJ_Depth", "SJ_Freq",
              "Ref_Motif", "Possivle_Alt_Motif","Possible_Alt_key", "Is_Canonical", "SJ_Overlap_Count", 
              "Mut_Pos", "Mut_Ref", "Mut_Alt", "Mut_Count", "Mut_Depth", "Mut_Freq",
              "Realign_No_SJ_Neg", "Realign_No_SJ_Pos", "Realign_Target_SJ_Neg", "Reaglin_Target_SJ_Pos",
              "Realign_Normal_SJ_Neg", "Realign_Normal_SJ_Pos"]
    print('\t'.join(header), file = hout)
    
    with open(input_file, 'r') as hin:
        next(hin)
        for line in hin:
            
            F = line.rstrip('\n').split('\t')
            if F[22] != "True": 
                print('\t'.join([F[0], F[1], F[2], F[6], F[7], F[8], F[9], F[10], F[11], F[12], F[13], F[14], F[15],
                             F[16], F[17], F[18], F[20], str(len(F[19])), F[21]]) +"\t-\t-\t-\t-\t-\t-", file = hout)                
            else:
                F = line.rstrip('\n').split('\t')
                mut_chr, junc_start, junc_end = F[0], int(F[1]), int(F[2])
                mut_pos, mut_ref, mut_alt = int(F[16]), F[17], F[18]
                if (F[6].endswith("5'SS") and F[7] == '+') or (F[6].endswith("3'SS") and F[7] == '-'):
                    junc_annotated = junc_end
                else:
                    junc_annotated = junc_start
    
                generate_template_seq(output_file + ".tmp.template.fa", ref_tb, junc_tb, mut_chr, mut_pos, mut_ref, mut_alt,
                    junc_start, junc_end, junc_annotated, template_size, genome_id, is_grc)
        
                extract_read_around_boundary(bam_file, output_file + ".tmp.read_seq.fa", mut_chr, mut_pos)
    
                type2count, ir_pos_sreads, target_pos_sreads = realign_main(output_file + ".tmp.read_seq.fa",
                    output_file + ".tmp.template.fa", 4 * template_size - score_margin)
    
                print('\t'.join([F[0], F[1], F[2], F[6], F[7], F[8], F[9], F[10], F[11], F[12], F[13], F[14], F[15],
                                 F[16], F[17], F[18], F[20], str(len(F[19])), F[21],
                                 str(type2count["no_splicing_negative"]), str(type2count["no_splicing_positive"]), 
                                 str(type2count["target_splicing_negative"]), str(type2count["target_splicing_positive"]),
                                 str(type2count["normal_splicing_negative"]), str(type2count["normal_splicing_positive"])]), file = hout)
                Path(output_file + ".tmp.template.fa").unlink()
                Path(output_file + ".tmp.read_seq.fa").unlink()

    ref_tb.close()
    junc_tb.close()
    
    Path(output_file + ".gencode.junc.bed.gz").unlink()
    Path(output_file + ".gencode.junc.bed.gz.tbi").unlink()

    hout.close()


if __name__ == "__main__":

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bam_file = sys.argv[3]
    reference = sys.argv[4]
    genome_id = sys.argv[5]
    is_grc = True if sys.argv[6] in ["True", "T", "TRUE"] else False

    #juncmut_realign(input_file, output_file, bam_file, reference, genome_id, is_grc, template_size = 10)
    juncmut_realign(input_file, output_file, bam_file, reference, genome_id, is_grc, template_size = 10, score_margin = 4)
    
"""
chr10      21659045        21661704        21659045        21661704        dummy   Intronic alternative 5'SS       -       108     184     0.5869565217391305      TCCCACCT        A.T.....        21661699:T:A,21661701:C:T       non-canonical        20      21661701        C       T       -     ,,....,,..,....TTtTt..  5       0.22727272727272727     True
chr5       124287120       124288581       124287118;124287120     124288579;124288581     dummy   Intronic alternative 3'SS       +       4       34      0.11764705882352941     ATTTCAG Y...AG. 124288580:C:A:124288581:A:G     CA|>AG*d     13      124288580       C       A       -       AAaA... 4       0.5714285714285714      True

""" 
