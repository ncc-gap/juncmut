#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
import random
import pysam
import annot_utils
import edlib

def change_base_check(tseq, tpos, tref, talt):
    # for debugging
    tseq_alt = tseq
    if tseq_alt[tpos] != tref:
        raise Exception("juncmut_realign.py: Inconsistent ref base")
    tseq_alt = tseq[:tpos] + talt + tseq[(tpos + 1):]
    return(tseq_alt)

def generate_template_seq(
        template_file, ref_tb, junc_tb, mut_chr, mut_pos, mut_ref, mut_alt, 
        junc_start, junc_end, junc_annotated, template_size):

    with open(template_file, 'w') as hout:
        # no splicing
        # without mutation
        ir_ref_seq = ref_tb.fetch(region = "%s:%d-%d" % (mut_chr, mut_pos - template_size, mut_pos + template_size - 1))
        print("No_SJ_Neg\t" + ir_ref_seq, file = hout)

        # with mutation
        ir_alt_seq = change_base_check(ir_ref_seq, template_size, mut_ref, mut_alt)
        print("No_SJ_Pos\t" + ir_alt_seq, file = hout)

        # target splicing
        # without mutation
        target_junc_ref_seq = ref_tb.fetch(region = "%s:%d-%d" % (mut_chr, (junc_start - 1) - (template_size - 1), junc_start - 1)) + \
            ref_tb.fetch(region = "%s:%d-%d" % (mut_chr, junc_end + 1, junc_end + 1 + template_size - 1))
        print("Target_SJ_Neg\t" + target_junc_ref_seq, file = hout)

        # with mutation
        if mut_pos >= (junc_start - 1) - (template_size - 1) and mut_pos <= junc_start - 1:
            target_junc_alt_seq = change_base_check(target_junc_ref_seq, mut_pos - ((junc_start - 1) - (template_size - 1)), mut_ref, mut_alt)
            print("Target_SJ_Pos\t" + target_junc_alt_seq, file = hout)
        elif mut_pos >= junc_end + 1 and mut_pos <= (junc_end + 1) + (template_size - 1):
            target_junc_alt_seq = change_base_check(target_junc_ref_seq, mut_pos - ((junc_end + 1) + (template_size - 1)) - 1, mut_ref, mut_alt)
            print("Target_SJ_Pos\t" + target_junc_alt_seq, file = hout)

        # normal splicing (multiple)
        records = junc_tb.fetch(region = "%s:%d-%d" % (mut_chr, mut_pos - template_size - 1, mut_pos + template_size))
        njind = 0
        for record_line in records:
            record = record_line.split('\t')
            rj_start, rj_end = int(record[1]) - 1, int(record[2])
            if (junc_annotated == junc_start and rj_start == junc_start and abs(rj_end - mut_pos) <= template_size) or \
                (junc_annotated == junc_end and rj_end == junc_end and abs(rj_start - 1 - mut_pos) <= template_size):
                # without mutation
                normal_junc_ref_seq = ref_tb.fetch(region = "%s:%d-%d" % (mut_chr, rj_start - 1 - (template_size - 1), rj_start - 1)) + \
                    ref_tb.fetch(region = "%s:%d-%d" % (mut_chr, rj_end + 1, rj_end + 1 + template_size - 1))
                print("Normal_SJ_Neg_%d\t%s" % (njind, normal_junc_ref_seq), file = hout)
                if mut_pos >= (rj_start - 1) - (template_size - 1) and mut_pos <= rj_start - 1:
                    normal_junc_alt_seq = change_base_check(normal_junc_ref_seq, mut_pos - (rj_start - 1 - (template_size - 1)), mut_ref, mut_alt)
                    print("Normal_SJ_Pos_%d\t%s" % (njind, normal_junc_alt_seq), file = hout)  
                elif mut_pos >= rj_end + 1 and mut_pos <= (rj_end + 1) + (template_size - 1):
                    normal_junc_alt_seq = change_base_check(normal_junc_alt_seq, mut_pos - (rj_end + 1 + template_size - 1) - 1, mut_ref, mut_alt)
                    print("Normal_SJ_Pos_%s\t%s" % (njind, normal_junc_alt_seq), file = hout)

            njind += 1

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

    return(check_flag)

def extract_read_around_boundary(bam_file, output_file, mut_chr, mut_pos, max_count = 10000):
    if bam_file.endswith('.bam'):
        bamfile = pysam.AlignmentFile(bam_file, 'rb')
    elif bam_file.endswith('.cram'):
        bamfile = pysam.AlignmentFile(bam_file, 'rc')
    
    read_count = bamfile.count(region = "%s:%d-%d" % (mut_chr, mut_pos, mut_pos), read_callback = check_read)
    read_ind_list = random.sample(range(read_count), min(max_count, read_count))
    read_inds = {}
    for i in read_ind_list:
        read_inds[i] = 1

    ind = -1
    with open(output_file, 'w') as hout:
        for read in bamfile.fetch(region = "%s:%d-%d" % (mut_chr, mut_pos, mut_pos)):
            if not check_read(read):
                continue

            ind += 1
            if not ind in read_inds:
                continue

            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]
            read_id = read.qname + '_1' if flags[6] == '1' else read.qname + '_2'
            print(">%s_%d_%d_%s_%d\n%s" % (read_id, read.reference_start + 1, read.reference_end, read.cigarstring, read.next_reference_start + 1, read.seq), file = hout)

    bamfile.close()

def process_cigar(res_cigar):
    #res_cigar = '16=1X1=1I1='
    # make a alignment object
    cigar_align = ''
    re_cigar = res_cigar.translate(str.maketrans({'=': ',=,','X': ',X,', 'I': ',I,', 'D': ',D,'}))
    cigar_list = re_cigar.split(',')
    del cigar_list[-1]
    for pnum in range(0, len(cigar_list), 2):
        cigar_align = cigar_align + str(str(cigar_list[pnum+1]) * int(cigar_list[pnum]))
    proc_cigar_list = list(cigar_align)
    
    for i in range(0, 5):
        pop_start = proc_cigar_list.pop(0) 
        pop_end = proc_cigar_list.pop(-1) 
        if pop_start == 'I' and pop_end == 'I':
            del proc_cigar_list[0]
            del proc_cigar_list[-1]
        elif pop_start == 'I' and pop_end != 'I':
            del proc_cigar_list[0]
        elif pop_start != 'I' and pop_end == 'I':
            del proc_cigar_list[-1]
        elif pop_start != 'I' and pop_end != 'I': continue
    
    mut_count = len([i for i in proc_cigar_list if i != '='])
    return(mut_count)

def juncmut_realign(input_file, output_file, bam_file, reference, genome_id, is_grc, template_size = 10):

    ref_tb = pysam.FastaFile(reference)

    if str(is_grc) == 'False':
        annot_utils.junction.make_junc_info(output_file + ".gencode.junc.bed.gz", "gencode", genome_id, False, False)
    else:
        annot_utils.junction.make_junc_info(output_file + ".gencode.junc.bed.gz", "gencode", genome_id, True, False)

    junc_tb = pysam.TabixFile(output_file + ".gencode.junc.bed.gz")
    sample_name = input_file.split("/")[-1].split('.')[0]

    # for each mutation in SJ.
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        output_header = ["Mut_key"] + csvreader.fieldnames + [
            "Sample", "Is_in_exon", "Realign_no_SJ_neg", "Realign_no_SJ_pos", "Realign_target_SJ_neg", "Reaglin_target_SJ_pos", "Realign_normal_SJ_neg", "Realign_normal_SJ_pos"
        ]
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames = output_header)
        csvwriter.writeheader()

        for csvobj in csvreader:
            output_obj = {}
            for key in csvobj:
                if key in output_header:
                    output_obj[key] = csvobj[key]

            (sj_key_chr, sj_key_start, sj_key_end) = csvobj["SJ_key"].split(",")
            sj_key_start = int(sj_key_start)
            sj_key_end = int(sj_key_end)
            output_obj["SJ_key"] = "%s:%d-%d" % (sj_key_chr, sj_key_start, sj_key_end)

            (mut_key_chr, mut_key_pos, mut_key_ref, mut_key_alt) = csvobj["Mut_key"].split(",")
            mut_key_pos = int(mut_key_pos)

            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]

            # [TODO] fix this item
            #output_obj["Mut_key"] = "%s,%s,%s,%s" % (csvobj["Chr"], csvobj["Mut_pos"], csvobj["Mut_ref"], csvobj["Mut_alt"])
            #output_obj["SJ_key"] = "%s:%s-%s" % (csvobj["Chr"], csvobj["Start"], csvobj["End"])
            output_obj["Sample"] = sample_name

            # key_readid_dict to remove read duplicates
            # Is a position of mutation in Exon or Intron
            is_exon = 'outside of motif'
            #o-->
            if splice_type == "Donor+":
                dis_mut_SJstart = mut_key_pos - sj_key_start
                if dis_mut_SJstart < 0 and dis_mut_SJstart > -3:
                    is_exon = 'exon'
                elif dis_mut_SJstart < 6 and dis_mut_SJstart >= 0:
                    is_exon = 'intron'
            #-->o
            elif splice_type == "Acceptor+":
                dis_mut_SJend = mut_key_pos - sj_key_end
                if dis_mut_SJend == 1:
                    is_exon = 'exon'
                elif dis_mut_SJend < 1 and dis_mut_SJend > -6:
                    is_exon = 'intron'
            #o<--       
            elif splice_type == "Acceptor-":
                dis_mut_SJstart = mut_key_pos - sj_key_start
                if dis_mut_SJstart == -1:
                    is_exon = 'exon'
                elif dis_mut_SJstart < 6 and dis_mut_SJstart >= 0:
                    is_exon = 'intron'
            #<--o
            elif splice_type == "Donor-": 
                dis_mut_SJend = mut_key_pos - sj_key_end
                if dis_mut_SJend > 0 and dis_mut_SJend < 3:
                    is_exon = 'exon'
                elif dis_mut_SJend < 1 and dis_mut_SJend > -6:
                    is_exon = 'intron'
            
            output_obj["Is_in_exon"] = is_exon

            #if RNA_mutation reads>=2 and Freq>=0.05, do realign.
            #if float(csvobj["Pileup_mut_count"])/float(csvobj["Pileup_depth"]) < mut_freq_thres or int(csvobj["Mut_count"]) < mut_num_thres:
            #    continue
            
            junc_annotated = sj_key_start
            if splice_type == "Donor+" or splice_type == "Acceptor-":
                junc_annotated = sj_key_end
            
            generate_template_seq(
                output_file + ".tmp.template.fa", ref_tb, junc_tb, mut_key_chr, mut_key_pos, mut_key_ref, mut_key_alt,
                sj_key_start, sj_key_end, junc_annotated, template_size)
            
            extract_read_around_boundary(bam_file, output_file + ".tmp.read_seq.fa", mut_key_chr, mut_key_pos)
            
            # query to dict
            d_query = {}
            with open(output_file + ".tmp.template.fa", 'r') as hin_template:
                for row in hin_template:
                    Q = row.rstrip('\n').split('\t')
                    d_query[Q[0]] = Q[1]

            # alignment
            d_realign = {'No_SJ_Neg': 0, 'No_SJ_Pos': 0, 'Target_SJ_Neg': 0, 'Target_SJ_Pos': 0, 'Normal_SJ_Neg': 0, 'Normal_SJ_Pos': 0}
            d_order = {'Normal_SJ_Neg': 0, 'No_SJ_Neg': 1, 'Target_SJ_Neg': 2, 'Normal_SJ_Pos': 3, 'No_SJ_Pos': 4, 'Target_SJ_Pos': 5}
            realign_result = "False"
            with open(output_file + ".tmp.read_seq.fa", 'r') as hin_read_seq:
                for row in hin_read_seq:
                    if row.startswith('>'):
                        # [TODO] Modify this section later
                        R = hin_read_seq.readline().rstrip('\n')
                        d_ed = {}
                        d_mut_count = {}
                        for key in d_query:
                            res = edlib.align(d_query[key], R, mode="HW", task="path")
                            res_cigar = res.get('cigar')
                            res_ed = res.get('editDistance')
                            # mut_count and edit distance
                            if res_ed <= 2:
                                # [TODO] rename mut_count
                                mut_count = process_cigar(res_cigar)
                                # [TODO] fix priority of if statement
                                # There may be multiple Normal SJ templates.  Choose template with smaller editDistance than 2.
                                if key.startswith('Normal_SJ_Neg_'):
                                    if 'Normal_SJ_Neg' in d_ed and d_ed['Normal_SJ_Neg'] <= res_ed: continue
                                    d_ed['Normal_SJ_Neg'] = res_ed
                                    
                                    if 'Normal_SJ_Neg' in d_mut_count and d_mut_count['Normal_SJ_Neg'] >= mut_count: continue
                                    d_mut_count['Normal_SJ_Neg'] = mut_count
                                
                                elif key.startswith('Normal_SJ_Pos_'):
                                    if 'Normal_SJ_Pos' in d_ed and d_ed['Normal_SJ_Pos'] <= res_ed: continue
                                    d_ed['Normal_SJ_Pos'] = res_ed 

                                    if 'Normal_SJ_Pos' in d_mut_count and d_mut_count['Normal_SJ_Pos'] >= mut_count: continue
                                    d_mut_count['Normal_SJ_Pos'] = mut_count
                                else:
                                    d_ed[key] = res_ed
                                    d_mut_count[key] = mut_count
                        
                        # choose minmum ed_distance and the key and count up.
                        # [TODO] add case, d_ed == {}: continue
                        min_k = [kv[0] for kv in d_ed.items() if kv[1] == min(d_ed.values())]
                        # [TODO] remove this case
                        if len(min_k) == 0: continue

                        min_k_sort = sorted(min_k, key=lambda x: d_order[x])
                        # if multiple mutation for the selected key is exist, no count-up.
                        if d_mut_count[min_k_sort[0]] == 0:
                            d_realign[min_k_sort[0]] += 1
                
                if d_realign['Target_SJ_Pos'] + d_realign['No_SJ_Pos']+ d_realign['Normal_SJ_Pos'] >=1:
                    realign_result = 'True'
            
            if realign_result == 'True':
                output_obj["Realign_no_SJ_neg"] = d_realign["No_SJ_Neg"]
                output_obj["Realign_no_SJ_pos"] = d_realign["No_SJ_Pos"]
                output_obj["Realign_target_SJ_neg"] = d_realign["Target_SJ_Neg"]
                output_obj["Reaglin_target_SJ_pos"] = d_realign["Target_SJ_Pos"]
                output_obj["Realign_normal_SJ_neg"] = d_realign["Normal_SJ_Neg"]
                output_obj["Realign_normal_SJ_pos"] = d_realign["Normal_SJ_Pos"]
                csvwriter.writerow(output_obj)

            os.remove(output_file + ".tmp.template.fa")
            os.remove(output_file + ".tmp.read_seq.fa")
            
    ref_tb.close()
    junc_tb.close()
    
    os.remove(output_file + ".gencode.junc.bed.gz")
    os.remove(output_file + ".gencode.junc.bed.gz.tbi")

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bam_file = sys.argv[3]
    reference = sys.argv[4]

    juncmut_realign(input_file, output_file, bam_file, reference, "hg38", "False", 10)
