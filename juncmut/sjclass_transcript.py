#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:01:08 2021

@author: Naoko Iida
gencode
wgEncodeGencodeBasicV39_hg38.txt.gz
CREATE TABLE `` (
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL, 0-start
  `cdsEnd` int(10) unsigned NOT NULL,   close
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `score` int(11) DEFAULT NULL,
  `name2` varchar(255) NOT NULL,
  `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `exonFrames` longblob NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`),
  KEY `name2` (`name2`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..

input_file
Mut_key = F[0] 
SJ_key = F[1]
    sj_chr = sj_elm[0]
    sj_pos_s = sj_elm[1].split('-')[0]
    sj_pos_e = sj_elm[1].split('-')[1]
Sample = F[2]
SJ_type = F[3]
SJ_strand = F[4]

list (tab separation)
[0]sample_name
[1]SJ.out.tab(full path)
[2]bam(full path)

"""

import json
import pysam
import csv

GENCODE_HEADER = [
  'bed_chrom',
  'bed_srat',
  'bed_end',
  'bin',
  'name',
  'chrom',
  'strand',
  'txStart',
  'txEnd',
  'cdsStart',
  'cdsEnd',
  'exonCount',
  'exonStarts',
  'exonEnds',
  'score',
  'name2',
  'cdsStartStat',
  'cdsEndStat',
  'exonFrames'
]

def longest_cds_tx(tx_info_list):
    if len(tx_info_list) == 0:
        return None

    model_tx_max_cds_length_list = []
    max_cds_length = 0
    for tx_obj in tx_info_list:
        gencode_exon_starts = list(map(int, tx_obj["Gencode_exon_starts"].rstrip(',').split(',')))
        gencode_exon_ends = list(map(int, tx_obj["Gencode_exon_ends"].rstrip(',').split(',')))
        
        cds_length = 0
        for n in range(0, len(gencode_exon_starts), 1):
            cds_length += gencode_exon_ends[n] - gencode_exon_starts[n]
        
        if max_cds_length < cds_length:
            max_cds_length = cds_length
            model_tx_max_cds_length_list = [tx_obj]
        elif max_cds_length == cds_length:
            model_tx_max_cds_length_list.append(tx_obj)

    if len(model_tx_max_cds_length_list) == 1:
        return model_tx_max_cds_length_list[0]

    # sort by Gencode_name
    model_tx_gencode_name_dict = {}
    for candidate in model_tx_max_cds_length_list:
        model_tx_gencode_name_dict[candidate['Gencode_name']] = candidate

    first_gencode_name = sorted(list(model_tx_gencode_name_dict.keys()))[0]
    return model_tx_gencode_name_dict[first_gencode_name]

def define_transcript(splice_type, juncmut_primary_sj_chr, juncmut_matching_ss, gencode_tb, mane_json):

    gencode_records = gencode_tb.fetch(region = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_matching_ss, juncmut_matching_ss + 1))
    if gencode_records == None:
        gencode_records = []

    tx_info_mane_list = []
    tx_info_mane_clinical_list = []
    tx_info_list = []
    for record in gencode_records:
        F = record.rstrip('\n').split('\t')
        gencode_chrom = F[GENCODE_HEADER.index("chrom")]
        # 0-base exon start == sj_pos
        gencode_exon_starts = list(map(int, F[GENCODE_HEADER.index("exonStarts")].rstrip(',').split(',')))
        # 0-base exon end == sj_pos +1
        gencode_exon_ends = list(map(int, F[GENCODE_HEADER.index("exonEnds")].rstrip(',').split(',')))
        gencode_cds_start = F[GENCODE_HEADER.index("cdsStart")]
        gencode_cds_end = F[GENCODE_HEADER.index("cdsEnd")]
        gencode_exon_count = int(F[GENCODE_HEADER.index("exonCount")])
        gencode_name2 = F[GENCODE_HEADER.index("name2")]
        gencode_name = F[GENCODE_HEADER.index("name")]
        gencode_exon_frames = F[GENCODE_HEADER.index("exonFrames")]

        juncmut_hijacked_ss = -1
        juncmut_hijacked_exon_num = -1
        if splice_type == "Donor+":
            if juncmut_matching_ss in gencode_exon_starts and juncmut_primary_sj_chr == gencode_chrom:
                index_start = gencode_exon_starts.index(juncmut_matching_ss)
                if index_start == 0:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_ends[index_start - 1] + 1
                    # be careful, the index_start is a 0-start count. juncmut_hijacked_exon_num is a 1-start count
                    juncmut_hijacked_exon_num = index_start - 1
        
        elif splice_type == "Donor-":
            if juncmut_matching_ss-1 in gencode_exon_ends and juncmut_primary_sj_chr == gencode_chrom:
                index_end = gencode_exon_ends.index(juncmut_matching_ss - 1)
                if index_end == len(gencode_exon_ends) - 1:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_starts[index_end + 1]
                    juncmut_hijacked_exon_num = index_end + 1
        
        elif splice_type == "Acceptor+":
            if juncmut_matching_ss-1 in gencode_exon_ends and juncmut_primary_sj_chr == gencode_chrom:
                index_end = gencode_exon_ends.index(juncmut_matching_ss - 1)
                if index_end == len(gencode_exon_ends) - 1:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_starts[index_end + 1]
                    juncmut_hijacked_exon_num = index_end + 1
        
        elif splice_type == "Acceptor-":
            if juncmut_matching_ss in gencode_exon_starts and juncmut_primary_sj_chr == gencode_chrom:
                index_start = gencode_exon_starts.index(juncmut_matching_ss)
                if index_start == 0:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_ends[index_start - 1] + 1
                    juncmut_hijacked_exon_num = index_start - 1
        
        if juncmut_hijacked_ss != -1 and juncmut_hijacked_exon_num != -1:
            tx_info = {
                "Juncmut_hijacked_SS": juncmut_hijacked_ss, 
                "Gencode_exon_frames": gencode_exon_frames, 
                "Gencode_name": gencode_name, 
                "Gencode_name2": gencode_name2, 
                "Juncmut_hijacked_exon_num": "%d/%d," % (juncmut_hijacked_exon_num, gencode_exon_count), 
                "Gencode_CDS_start": gencode_cds_start, 
                "Gencode_CDS_end": gencode_cds_end, 
                "Gencode_exon_starts": F[GENCODE_HEADER.index("exonStarts")], 
                "Gencode_exon_ends": F[GENCODE_HEADER.index("exonEnds")],
                "MANE": "-"
            }
            if gencode_name in mane_json:
                if mane_json[gencode_name]["tag"] == "MANE_Select":
                    tx_info["MANE"] = "MANE_Select"
                    tx_info_mane_list.append(tx_info)
                elif mane_json[gencode_name]["tag"] == "MANE_Plus_Clinical":
                    tx_info["MANE"] = "MANE_Plus_Clinical"
                    tx_info_mane_clinical_list.append(tx_info)
                else:
                    tx_info["MANE"] = "not_MANE_tag"
                    tx_info_list.append(tx_info)
            else:
                tx_info["MANE"] = "unmatched_MANE"
                tx_info_list.append(tx_info)

    model_tx = longest_cds_tx(tx_info_mane_list)
    if model_tx is None:
        model_tx = longest_cds_tx(tx_info_mane_clinical_list)
    if model_tx is None:
        model_tx = longest_cds_tx(tx_info_list)

    return model_tx

def sjclass_transcript(input_file, output_file, gencode, mane):

    gencode_tb = pysam.TabixFile(gencode)

    with open(mane, 'r') as hin_mane:
        mane_json = json.load(hin_mane)

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=[
            "Mut_key", "SJ_key", "Sample", "SJ_type", "SJ_strand", "Is_GTAG_creation", "Is_in_new_exon",
            "Realign_no_SJ_neg", "Realign_no_SJ_pos", "Realign_target_SJ_neg", "Reaglin_target_SJ_pos", "Realign_normal_SJ_neg", "Realign_normal_SJ_pos",
            "gnomAD", "gnomAD_AF", 
            "Gencode_name", "Gencode_name2", "Gencode_CDS_start", "Gencode_CDS_end", "Gencode_exon_starts", "Gencode_exon_ends", "Gencode_exon_frames", "Gene_type",
            "Juncmut_hijacked_exon_num","Juncmut_hijacked_SS", "Juncmut_primary_SS", "Juncmut_matching_SS", "SS_shift1", "SS_shift2", "MANE"
        ])
        csvwriter.writeheader()

        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            (juncmut_primary_sj_chr, juncmut_primary_sj_pos) = csvobj["SJ_key"].split(":")
            (juncmut_primary_sj_start, juncmut_primary_sj_end) = list(map(int, juncmut_primary_sj_pos.split('-')))
            juncmut_gene_strand = csvobj["SJ_strand"]
            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]
            sample = csvobj["Sample"]

            # o--->
            if splice_type == "Donor+":
                juncmut_matching_ss = juncmut_primary_sj_end
                juncmut_primary_ss = juncmut_primary_sj_start
            # <---o
            elif splice_type == "Donor-":
                juncmut_matching_ss = juncmut_primary_sj_start
                juncmut_primary_ss = juncmut_primary_sj_end
            # --->o
            elif splice_type == "Acceptor+":
                juncmut_matching_ss = juncmut_primary_sj_start
                juncmut_primary_ss = juncmut_primary_sj_end
            # o<---
            elif splice_type == "Acceptor-":
                juncmut_matching_ss = juncmut_primary_sj_end
                juncmut_primary_ss = juncmut_primary_sj_start

            model_tx = define_transcript(splice_type, juncmut_primary_sj_chr, juncmut_matching_ss, gencode_tb, mane_json)
            if model_tx is None: continue
            
            # The ss_shift1 from the corresponding normal splice site. Minus indicates 5' side on the genome. Plus indicates 3'SS.
            ss_shift1 = juncmut_primary_ss - model_tx["Juncmut_hijacked_SS"]
            ss_shift2 = ss_shift1
            # The ss_shift1 from the corresponding normal splice site. Minus indicates the exon side on the pre-mRNA. Minus indicates exon side
            if splice_type in ["Donor-", "Acceptor+"]:
                ss_shift2 = -ss_shift1

            # Is_GTAG_creation
            Is_GTAG_creation = "True"
            if csvobj["Is_GT_AG"].startswith("non-"):
                Is_GTAG_creation = "False"

            # gene_type
            gene_type = "Coding"
            if model_tx["Gencode_CDS_start"] == model_tx["Gencode_CDS_end"]:
                gene_type = "Non-coding"

            new_csvobj = {
                "Mut_key": csvobj["Mut_key"],
                "SJ_key": csvobj["SJ_key"],
                "Sample": csvobj["Sample"],
                "SJ_type": csvobj["Created_motif"],
                "SJ_strand": csvobj["SJ_strand"],
                "Is_GTAG_creation": Is_GTAG_creation,
                "Is_in_new_exon": csvobj["Is_in_exon"],
                "Realign_no_SJ_neg": csvobj["Realign_no_SJ_neg"],
                "Realign_no_SJ_pos": csvobj["Realign_no_SJ_pos"],
                "Realign_target_SJ_neg": csvobj["Realign_target_SJ_neg"],
                "Reaglin_target_SJ_pos": csvobj["Reaglin_target_SJ_pos"],
                "Realign_normal_SJ_neg": csvobj["Realign_normal_SJ_neg"],
                "Realign_normal_SJ_pos": csvobj["Realign_normal_SJ_pos"],
                "gnomAD": csvobj["gnomAD"],
                "gnomAD_AF": csvobj["gnomAD_AF"],
                "Gencode_name": model_tx["Gencode_name"],
                "Gencode_name2": model_tx["Gencode_name2"],
                "Gencode_CDS_start": model_tx["Gencode_CDS_start"],
                "Gencode_CDS_end": model_tx["Gencode_CDS_end"],
                "Gencode_exon_starts": model_tx["Gencode_exon_starts"],
                "Gencode_exon_ends": model_tx["Gencode_exon_ends"],
                "Gencode_exon_frames": model_tx["Gencode_exon_frames"],
                "Gene_type": gene_type,
                "Juncmut_hijacked_exon_num": model_tx["Juncmut_hijacked_exon_num"],
                "Juncmut_hijacked_SS": model_tx["Juncmut_hijacked_SS"],
                "Juncmut_primary_SS": juncmut_primary_ss,
                "Juncmut_matching_SS": juncmut_matching_ss,
                "SS_shift1": ss_shift1,
                "SS_shift2": ss_shift2,
                "MANE": model_tx["MANE"],
            }
            csvwriter.writerow(new_csvobj)

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gencode = sys.argv[3]
    mane = sys.argv[4]

    sjclass_transcript(input_file, output_file, gencode, mane)
