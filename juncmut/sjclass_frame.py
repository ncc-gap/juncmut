#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:01:08 2021

@author: Naoko Iida
fetch 0-start open-close type(exclude-include).
juncmut_primary_ss : junction position of intron side.

dist = mut - norm (distance to mut from normal position.)

txStart
txEnd
cdsStart
cdsEnd
"""

import pysam
import csv
import juncmut.protein

DICT_CODON = {"TTT":"F", "TTC":"F", "TCT":"S", "TCC":"S", "TAT":"Y", "TAC":"Y", "TGT":"C", "TGC":"C",
          "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*", "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W", 
          "CTT":"L", "CTC":"L", "CCT":"P", "CCC":"P", "CAT":"H", "CAC":"H", "CGT":"R", "CGC":"R", 
          "CTA":"L", "CTG":"L", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGG":"R", 
          "ATT":"I", "ATC":"I", "ACT":"T", "ACC":"T", "AAT":"N", "AAC":"N", "AGT":"S", "AGC":"S", 
          "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R", "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R", 
          "GTT":"V", "GTC":"V", "GCT":"A", "GCC":"A", "GAT":"D", "GAC":"D", "GGT":"G", "GGC":"G", 
          "GTA":"V", "GTG":"V", "GCA":"A", "GCG":"A", "GAA":"E", "GAG":"E", "GGA":"G", "GGG":"G"}

def translate(org_seq, is_inframe, is_ptc, is_reference, ref_cds_exon_len, ref_cds_len, last_exon_len):

    aa_seq = ""

    is_inframe_ptc = not is_reference and is_inframe == "In-frame" and is_ptc == "PTC"
    reference_stop_codon_pos = 0
    if is_inframe_ptc:
        reference_stop_codon_pos = ref_cds_len + len(org_seq) - ref_cds_exon_len

    for seq_index in range(0, len(org_seq), 3):
        codon = org_seq[seq_index:seq_index+3]
        aa = DICT_CODON.get(codon, '-')
        aa_seq += aa
        if aa == "*":
            if not is_inframe_ptc:
                break
        if is_inframe_ptc and reference_stop_codon_pos == seq_index + 3:
            break

    return aa_seq

def inspection(q_target_seq, ref_cds_exon_len, ref_cds_len, last_exon_len):

    aa_seq = ""
    is_ptc = "NA"
    is_nmd = "NA"
    non_nmd_region = last_exon_len + 50
    # reference stop codon
    reference_stop_codon_pos = ref_cds_len + len(q_target_seq) - ref_cds_exon_len

    #remove the actual last codon
    for seq_index in range(0, len(q_target_seq), 3):
        codon = str(q_target_seq[seq_index:seq_index+3])
        aa = DICT_CODON.get(codon, '-')
        aa_seq += aa
        if aa == "*":
            # 0aa-start +3 for adjust
            if reference_stop_codon_pos == seq_index + 3:
                is_ptc = "Survived"
            else:
                is_ptc = "PTC" 
                if len(q_target_seq) - seq_index > non_nmd_region:
                    is_nmd = "TRUE"
                else:
                    is_nmd = "FALSE"
            break

    if is_ptc == "NA":
        is_ptc = "Disappeared"
    return(is_ptc, is_nmd)

def revcomp(seq):
    trans = str.maketrans('ATGCatgc', 'TACGTAGC')
    rev = seq[::-1]
    revcomp = rev.translate(trans)
    return(revcomp)

def frame_shift(length):
    q, mod = divmod(length, 3)
    if mod == 0:
        frame = "In-frame"
    else: frame = "Frameshift"
    
    return(frame)


def convert_aa (seq, is_inframe, is_ptc):
    import re
    seq_dst = re.sub("^p\\.\\(", "", seq)
    seq_dst = re.sub("\\)$", "", seq_dst)
   
    seq_dst = re.sub("Ala", "A", seq_dst)
    seq_dst = re.sub("Arg", "R", seq_dst)
    seq_dst = re.sub("Asn", "N", seq_dst)
    seq_dst = re.sub("Asp", "D", seq_dst)
    seq_dst = re.sub("Cys", "C", seq_dst)
    seq_dst = re.sub("Gln", "Q", seq_dst)
    seq_dst = re.sub("Glu", "E", seq_dst)
    seq_dst = re.sub("Gly", "G", seq_dst)
    seq_dst = re.sub("His", "H", seq_dst)
    seq_dst = re.sub("Ile", "I", seq_dst)
    seq_dst = re.sub("Leu", "L", seq_dst)
    seq_dst = re.sub("Lys", "K", seq_dst)
    seq_dst = re.sub("Met", "M", seq_dst)
    seq_dst = re.sub("Phe", "F", seq_dst)
    seq_dst = re.sub("Pro", "P", seq_dst)
    seq_dst = re.sub("Ser", "S", seq_dst)
    seq_dst = re.sub("Thr", "T", seq_dst)
    seq_dst = re.sub("Trp", "W", seq_dst)
    seq_dst = re.sub("Tyr", "Y", seq_dst)
    seq_dst = re.sub("Val", "V", seq_dst)
    seq_dst = re.sub("Ter", "*", seq_dst)
    
    result_seq = ""
    p = re.compile(r"ins[A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|S|T|W|Y|V]{3,}")
    for item in seq_dst.replace("ins", ";ins").split(";"):
        if re.search(p, item):
            match = re.search(p, item).group()
            length_of_insertion = len(match) - 3
            replacement = f"insX[{length_of_insertion}]"
            result_seq += re.sub(p, replacement, item)
        else:
            result_seq += item

    if is_inframe == "In-frame" and is_ptc == "PTC":
        result_seq = result_seq.split("*")[0] + "*"
    
    return result_seq

def sjclass_frame(input_file, output_file, reference):
    
    ref = pysam.FastaFile(reference) # 0-coordinate
          
    # for each row.
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Mut_pos_in_motif", "Is_coding_SJ", "Start_codon_disruption_num", "Stop_codon_disruption_num", "Is_inframe", "Is_PTC", "Is_NMD",
            "Protein_description_change", "Protein_description_change_short", "Protein_description_first_pos", "Protein_description_last_pos_ref", "Protein_description_last_pos_target", "Protein_description_change_mut", "Protein_description_change_mut_short"
        ])
        csvwriter.writeheader()
        for csvobj in csvreader:
            ## set factors ##################################################
            mut_pos_in_motif = "NA"
            is_coding_sj = "Coding"
            is_inframe = "NA"
            is_ptc = "NA"
            is_nmd = "NA"
            protein_description_change = "NA"
            protein_description_change_short = "NA"
            protein_description_first_pos = "NA"
            protein_description_last_pos_ref = "NA"
            protein_description_last_pos_target = "NA"
            protein_description_change_mut = "NA"
            protein_description_change_mut_short = "NA"
            start_codon_disruption_num = 0
            stop_codon_disruption_num = 0

            csvobj["Mut_pos_in_motif"] = "NA"
            csvobj["Is_coding_SJ"] = "NA"
            csvobj["Start_codon_disruption_num"] = "NA"
            csvobj["Stop_codon_disruption_num"] = "NA"
            csvobj["Is_inframe"] = is_inframe
            csvobj["Is_PTC"] = is_ptc
            csvobj["Is_NMD"] = is_nmd
            csvobj["Protein_description_change"] = protein_description_change
            csvobj["Protein_description_change_short"] = protein_description_change_short
            csvobj["Protein_description_first_pos"] = protein_description_first_pos
            csvobj["Protein_description_last_pos_ref"] = protein_description_last_pos_ref
            csvobj["Protein_description_last_pos_target"] = protein_description_last_pos_target
            csvobj["Protein_description_change_mut"] = protein_description_change_mut
            csvobj["Protein_description_change_mut_short"] = protein_description_change_mut_short
            if csvobj["Gencode_exon_starts"] == "NA":
                csvwriter.writerow(csvobj)
                continue

            mut_key_pos = int(csvobj["Mut_key"].split(',')[1])
            mut_key_allele = csvobj["Mut_key"].split(',')[3]
            (sj_key_chr, sj_key_pos) = csvobj["SJ_key"].split(":")
            (sj_key_pos_start, sj_key_pos_end) = list(map(int, sj_key_pos.split('-')))

            juncmut_predicted_splicing_type = csvobj["Juncmut_predicted_splicing_type"]
            juncmut_primary_ss = int(csvobj["Juncmut_primary_SS"])
            
            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]
            sj_strand = csvobj["SJ_strand"]

            gencode_exon_starts = list(map(int, csvobj["Gencode_exon_starts"].rstrip(',').split(',')))
            gencode_exon_ends = list(map(int, csvobj["Gencode_exon_ends"].rstrip(',').split(',')))
            gencode_cds_start = int(csvobj["Gencode_CDS_start"].rstrip(','))
            gencode_cds_end = int(csvobj["Gencode_CDS_end"].rstrip(','))

            ## mut position in the motif ###################################################
            if splice_type in ["Donor+", "Acceptor-"]:
                mut_pos_in_motif = mut_key_pos - juncmut_primary_ss
            elif splice_type in ["Donor-", "Acceptor+"]:
                mut_pos_in_motif = juncmut_primary_ss - mut_key_pos
            else:
                mut_pos_in_motif = "NA"

            ## the gene region of sj

            if gencode_cds_start == gencode_cds_end:
                is_coding_sj = "Non-coding"
            else:
                if sj_strand == "+":
                    if sj_key_pos_end <= gencode_cds_start:
                        is_coding_sj = "5'UTR"
                    elif sj_key_pos_start > gencode_cds_end:
                        is_coding_sj = "3'UTR"
                    else:
                        if sj_key_pos_start-1 <= gencode_cds_start < sj_key_pos_end:
                            start_codon_disruption_num = sj_key_pos_end - gencode_cds_start
                        if sj_key_pos_start <= gencode_cds_end <= sj_key_pos_end:
                            stop_codon_disruption_num = gencode_cds_end + 1 - sj_key_pos_start

                        if 3 <= start_codon_disruption_num: 
                            is_coding_sj = "Start_codon_disruption1"
                        elif 3 <= stop_codon_disruption_num:
                            is_coding_sj = "Stop_codon_disruption1"

                else:
                    if sj_key_pos_start > gencode_cds_end:
                        is_coding_sj = "5'UTR"
                    elif sj_key_pos_end <= gencode_cds_start:
                        is_coding_sj = "3'UTR"
                    else:
                        if sj_key_pos_start <= gencode_cds_end <= sj_key_pos_end:
                            start_codon_disruption_num = gencode_cds_end + 1 - sj_key_pos_start
                        if sj_key_pos_start-1 <= gencode_cds_start < sj_key_pos_end:
                            stop_codon_disruption_num = sj_key_pos_end - gencode_cds_start

                        if 3 <= start_codon_disruption_num:
                            is_coding_sj = "Start_codon_disruption1"
                        elif 3 <= stop_codon_disruption_num:
                            is_coding_sj = "Stop_codon_disruption1"


            if is_coding_sj == "Coding" \
            and juncmut_predicted_splicing_type not in ["-", "Ambiguous_outgene", "Ambiguous_termination"]:

                (hijacked_exon_num, gencode_exon_count) = list(map(int, csvobj["Juncmut_hijacked_exon_num"].rstrip(',').split('/')))
                (closed_exon_num, gencode_exon_count) = list(map(int, csvobj["Closed_exon_num"].rstrip(',').split('/')))
                num_skipped_exon = int(csvobj["Num_skipped_exon"])
                juncmut_secondary_ss = csvobj["Juncmut_secondary_SS"]
                is_secondary_ss = False
                if juncmut_secondary_ss != "NA":
                    is_secondary_ss = True
                    juncmut_secondary_ss = int(juncmut_secondary_ss)

                # check the exon with gencode_cds_start and gencode_cds_end
                for exon_num in range(0, gencode_exon_count, 1):
                    if gencode_exon_starts[exon_num] <= gencode_cds_start <= gencode_exon_ends[exon_num]:
                        gencode_cds_start_exon_num = exon_num
                    if gencode_exon_starts[exon_num] <= gencode_cds_end <= gencode_exon_ends[exon_num]:
                        gencode_cds_end_exon_num = exon_num

                # new_exon_pos_list . start codon to exon end. (sj_strand dependent)
                new_exon_starts =[]
                new_exon_ends =[]
                new_exon_shift_exon_num = 0

                if sj_strand == "+":
                    new_exon_shift_exon_num = gencode_cds_start_exon_num
                    for exon_num in range(gencode_cds_start_exon_num, gencode_exon_count, 1):
                        if exon_num == gencode_cds_start_exon_num:
                            new_exon_starts.append(gencode_cds_start)
                            new_exon_ends.append(gencode_exon_ends[exon_num])
                        elif exon_num == gencode_cds_end_exon_num:
                            new_exon_starts.append(gencode_exon_starts[exon_num])
                            new_exon_ends.append(gencode_exon_ends[exon_num])
                        else:
                            new_exon_starts.append(gencode_exon_starts[exon_num])
                            new_exon_ends.append(gencode_exon_ends[exon_num])
                else:
                    for exon_num in range(0, gencode_cds_end_exon_num+1, 1):
                        if exon_num == gencode_cds_start_exon_num:
                            new_exon_starts.append(gencode_exon_starts[exon_num])
                            new_exon_ends.append(gencode_exon_ends[exon_num])
                        elif exon_num == gencode_cds_end_exon_num:
                            new_exon_starts.append(gencode_exon_starts[exon_num])
                            new_exon_ends.append(gencode_cds_end)
                        else:
                            new_exon_starts.append(gencode_exon_starts[exon_num])
                            new_exon_ends.append(gencode_exon_ends[exon_num])

                # calculate reference cds+3'utr-exon length
                ref_cds_exon_len = 0
                tmp_ref_cds_exon = ""
                tmp_ref_cds_exon_mut = ""
                is_inside_exon = False
                for exon_num in range(0, len(new_exon_starts), 1):
                    ref_cds_exon_len += new_exon_ends[exon_num] - new_exon_starts[exon_num]
                    read = ref.fetch(sj_key_chr,new_exon_starts[exon_num], new_exon_ends[exon_num])
                    tmp_ref_cds_exon += read

                    if new_exon_starts[exon_num] < mut_key_pos <= new_exon_ends[exon_num]:
                        index_mut_key_pos = mut_key_pos - new_exon_starts[exon_num] - 1
                        tmp_ref_cds_exon_mut += read[:index_mut_key_pos] + mut_key_allele + read[index_mut_key_pos + 1:]
                        is_inside_exon = True

                        mut_key_ref = csvobj["Mut_key"].split(',')[2]
                        if read[index_mut_key_pos] != mut_key_ref:
                            raise Exception("%s is not match mut_key_ref: %d" % (read[index_mut_key_pos], mut_key_ref))

                    else:
                        tmp_ref_cds_exon_mut += read

                if ref_cds_exon_len != len(tmp_ref_cds_exon):
                    raise Exception("length tmp_ref_cds_exon: %d is not match ref_cds_exon_len: %d" % (len(tmp_ref_cds_exon), ref_cds_exon_len))

                if sj_strand == "+":
                    ref_cds_exon = tmp_ref_cds_exon.upper()
                    ref_cds_exon_mut = tmp_ref_cds_exon_mut.upper()
                else:
                    ref_cds_exon = revcomp(tmp_ref_cds_exon.upper())
                    ref_cds_exon_mut = revcomp(tmp_ref_cds_exon_mut.upper())

                # calculate reference cds length
                ref_cds_len = 0
                if gencode_cds_start_exon_num == gencode_cds_end_exon_num:
                    ref_cds_len = gencode_cds_end - gencode_cds_start
                else:
                    for exon_num in range(gencode_cds_start_exon_num, gencode_cds_end_exon_num + 1, 1):
                        if exon_num == gencode_cds_start_exon_num:
                            ref_cds_len = ref_cds_len + gencode_exon_ends[exon_num] - gencode_cds_start
                        elif exon_num == gencode_cds_end_exon_num:
                            ref_cds_len = ref_cds_len + gencode_cds_end - gencode_exon_starts[exon_num]
                        else:
                            ref_cds_len = ref_cds_len + gencode_exon_ends[exon_num] - gencode_exon_starts[exon_num]

                target_seq = ""
                first_exon_seq = ""
                last_exon_seq = ""

                ## Main process ####
                if splice_type in ["Donor+", "Acceptor-"]:

                    for exon_num in range(0, len(new_exon_starts)):

                        if num_skipped_exon > 0 and closed_exon_num < exon_num + new_exon_shift_exon_num <= hijacked_exon_num:
                            continue

                        get_seq = ""
                        if exon_num + new_exon_shift_exon_num == closed_exon_num:
                            if juncmut_predicted_splicing_type == "Cryptic_exon":
                                if juncmut_secondary_ss < mut_key_pos < juncmut_primary_ss:
                                    intron_seq = ref.fetch(sj_key_chr, juncmut_secondary_ss, mut_key_pos-1) \
                                        + mut_key_allele \
                                        + ref.fetch(sj_key_chr,mut_key_pos, juncmut_primary_ss-1)
                                else:
                                    intron_seq = ref.fetch(sj_key_chr, juncmut_secondary_ss, juncmut_primary_ss-1)

                                exon_seq = ref.fetch(sj_key_chr,new_exon_starts[exon_num], new_exon_ends[exon_num])
                                get_seq = exon_seq + intron_seq
                                if exon_num == 0:
                                    first_exon_seq = exon_seq

                            elif juncmut_predicted_splicing_type == "Lengthening_exon" or juncmut_predicted_splicing_type == "Shortening_exon":
                                if mut_key_pos < juncmut_primary_ss:
                                    if new_exon_starts[exon_num] < mut_key_pos-1:
                                        get_seq = ref.fetch(sj_key_chr, new_exon_starts[exon_num], mut_key_pos-1)
                                    get_seq +=  mut_key_allele + ref.fetch(sj_key_chr, mut_key_pos, juncmut_primary_ss-1)
                                else:
                                    get_seq = ref.fetch(sj_key_chr, new_exon_starts[exon_num], juncmut_primary_ss-1)

                        else:
                            get_seq = ref.fetch(sj_key_chr,new_exon_starts[exon_num], new_exon_ends[exon_num])

                        if first_exon_seq == "":
                            first_exon_seq = get_seq
                        target_seq += get_seq

                    last_exon_seq = get_seq

                # 5'SS - <------o
                elif splice_type in ["Donor-", "Acceptor+"]:

                    for exon_num in range(0, len(new_exon_starts)):

                        if num_skipped_exon > 0 and hijacked_exon_num <= exon_num + new_exon_shift_exon_num < closed_exon_num:
                            continue

                        if exon_num + new_exon_shift_exon_num == closed_exon_num:
                            if juncmut_predicted_splicing_type == "Cryptic_exon":
                                if juncmut_primary_ss < mut_key_pos < juncmut_secondary_ss:
                                    intron_seq = ref.fetch(sj_key_chr, juncmut_primary_ss, mut_key_pos-1) \
                                        + mut_key_allele \
                                        + ref.fetch(sj_key_chr,mut_key_pos, juncmut_secondary_ss-1)
                                else:
                                    intron_seq = ref.fetch(sj_key_chr, juncmut_primary_ss, juncmut_secondary_ss-1)

                                exon_seq = ref.fetch(sj_key_chr,new_exon_starts[exon_num], new_exon_ends[exon_num])
                                get_seq = intron_seq + exon_seq
                                if exon_num == len(new_exon_starts) - 1:
                                    last_exon_seq = exon_seq

                            elif juncmut_predicted_splicing_type == "Lengthening_exon" or juncmut_predicted_splicing_type == "Shortening_exon":
                                if juncmut_primary_ss < mut_key_pos:
                                    get_seq = ref.fetch(sj_key_chr, juncmut_primary_ss, mut_key_pos-1) + mut_key_allele
                                    if mut_key_pos < new_exon_ends[exon_num]:
                                        get_seq += ref.fetch(sj_key_chr,mut_key_pos, new_exon_ends[exon_num])

                                else:
                                    get_seq = ref.fetch(sj_key_chr, juncmut_primary_ss, new_exon_ends[exon_num])

                        else:
                            get_seq = ref.fetch(sj_key_chr,new_exon_starts[exon_num], new_exon_ends[exon_num])

                        if exon_num == 0:
                            first_exon_seq = get_seq
                        target_seq += get_seq

                    if last_exon_seq == "":
                        last_exon_seq = get_seq

                if sj_strand == "+":
                    q_target_seq = target_seq.upper()
                    last_exon_len = len(last_exon_seq)
                else:
                    target_seq_tmp = target_seq.upper()
                    q_target_seq = revcomp(target_seq_tmp)
                    last_exon_len = len(first_exon_seq)

                is_ptc, is_nmd = inspection(q_target_seq, ref_cds_exon_len, ref_cds_len, last_exon_len)
                is_inframe = frame_shift(len(target_seq) - ref_cds_exon_len)

                if sj_strand == "+":
                    if start_codon_disruption_num in [1, 2] and aa_seq_ignore_stop[0] != "M":
                        is_coding_sj = "Start_codon_disruption2"
                    elif stop_codon_disruption_num in [1, 2] and aa_seq_ignore_stop[-1] != "*":
                        is_coding_sj = "Stop_codon_disruption2"

                else:
                    if start_codon_disruption_num in [1, 2] and aa_seq_ignore_stop[0] != "M":
                        is_coding_sj = "Start_codon_disruption2"
                    elif stop_codon_disruption_num in [1, 2] and aa_seq_ignore_stop[-1] != "*":
                        is_coding_sj = "Stop_codon_disruption2"

                q_target_seq_translate = translate(q_target_seq, is_inframe, is_ptc, False, ref_cds_exon_len, ref_cds_len, last_exon_len)
                ref_cds_exon_translate = translate(ref_cds_exon, None, None, True, 0, 0, 0)
                #result1 = protein_description(ref_cds_len, ref_cds_exon_translate, q_target_seq_translate)
                if is_inframe == "In-frame": 
                    result1 = juncmut.protein.in_frame_description(ref_cds_exon_translate, q_target_seq_translate)
                else:
                    result1 = juncmut.protein.out_of_frame_description(ref_cds_exon_translate, q_target_seq_translate)

                (protein_description_change, protein_description_first_pos, protein_description_last_pos_ref, protein_description_last_pos_target) = result1
                protein_description_change_short = convert_aa(protein_description_change, is_inframe, is_ptc)

                if is_inside_exon:
                    ref_cds_exon_mut_translate = translate(ref_cds_exon_mut, is_inframe, is_ptc, False, ref_cds_exon_len, ref_cds_len, last_exon_len)
                    #result2 = protein_description(ref_cds_len, ref_cds_exon_translate, ref_cds_exon_mut_translate)
                    result2 = juncmut.protein.in_frame_description(ref_cds_exon_translate, ref_cds_exon_mut_translate)
                    
                    protein_description_change_mut = result2[0]
                    protein_description_change_mut_short = convert_aa(protein_description_change_mut, is_inframe, is_ptc)

            csvobj["Mut_pos_in_motif"] = mut_pos_in_motif
            csvobj["Is_coding_SJ"] = is_coding_sj
            csvobj["Start_codon_disruption_num"] = start_codon_disruption_num
            csvobj["Stop_codon_disruption_num"] = stop_codon_disruption_num
            csvobj["Is_inframe"] = is_inframe
            csvobj["Is_PTC"] = is_ptc
            csvobj["Is_NMD"] = is_nmd
            csvobj["Protein_description_change"] = protein_description_change
            csvobj["Protein_description_change_short"] = protein_description_change_short
            csvobj["Protein_description_first_pos"] = protein_description_first_pos
            csvobj["Protein_description_last_pos_ref"] = protein_description_last_pos_ref
            csvobj["Protein_description_last_pos_target"] = protein_description_last_pos_target
            csvobj["Protein_description_change_mut"] = protein_description_change_mut
            csvobj["Protein_description_change_mut_short"] = protein_description_change_mut_short

            csvwriter.writerow(csvobj)

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reference = sys.argv[3]

    sjclass_frame(input_file, output_file, reference)
