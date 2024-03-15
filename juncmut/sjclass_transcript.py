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

def define_transcript(splice_type, juncmut_primary_sj_chr, juncmut_matching_ss, gencode_tb, transcript):

    gencode_records = gencode_tb.fetch(region = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_matching_ss, juncmut_matching_ss + 1))
    if gencode_records == None:
        return None

    tx_info = None
    for record in gencode_records:
        F = record.rstrip('\n').split('\t')
        gencode_chrom = F[GENCODE_HEADER.index("chrom")]
        if juncmut_primary_sj_chr != gencode_chrom:
            continue

        if transcript != F[GENCODE_HEADER.index("name")]:
            continue

        # 0-base exon start == sj_pos
        gencode_exon_starts = list(map(int, F[GENCODE_HEADER.index("exonStarts")].rstrip(',').split(',')))
        # 0-base exon end == sj_pos +1
        gencode_exon_ends = list(map(int, F[GENCODE_HEADER.index("exonEnds")].rstrip(',').split(',')))
        gencode_cds_start = F[GENCODE_HEADER.index("cdsStart")]
        gencode_cds_end = F[GENCODE_HEADER.index("cdsEnd")]
        gencode_exon_count = int(F[GENCODE_HEADER.index("exonCount")])
        gencode_exon_frames = F[GENCODE_HEADER.index("exonFrames")]

        juncmut_hijacked_ss = -1
        juncmut_hijacked_exon_num = -1
        if splice_type == "Donor+":
            if juncmut_matching_ss in gencode_exon_starts:
                index_start = gencode_exon_starts.index(juncmut_matching_ss)
                if index_start == 0:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_ends[index_start - 1] + 1
                    # be careful, the index_start is a 0-start count. juncmut_hijacked_exon_num is a 1-start count
                    juncmut_hijacked_exon_num = index_start - 1
        
        elif splice_type == "Donor-":
            if juncmut_matching_ss-1 in gencode_exon_ends:
                index_end = gencode_exon_ends.index(juncmut_matching_ss - 1)
                if index_end == len(gencode_exon_ends) - 1:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_starts[index_end + 1]
                    juncmut_hijacked_exon_num = index_end + 1
        
        elif splice_type == "Acceptor+":
            if juncmut_matching_ss-1 in gencode_exon_ends:
                index_end = gencode_exon_ends.index(juncmut_matching_ss - 1)
                if index_end == len(gencode_exon_ends) - 1:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_starts[index_end + 1]
                    juncmut_hijacked_exon_num = index_end + 1
        
        elif splice_type == "Acceptor-":
            if juncmut_matching_ss in gencode_exon_starts:
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
                "Juncmut_hijacked_exon_num": "%d/%d," % (juncmut_hijacked_exon_num, gencode_exon_count), 
                "Gencode_CDS_start": gencode_cds_start, 
                "Gencode_CDS_end": gencode_cds_end, 
                "Gencode_exon_starts": F[GENCODE_HEADER.index("exonStarts")], 
                "Gencode_exon_ends": F[GENCODE_HEADER.index("exonEnds")],
            }

    return tx_info

def sjclass_transcript(input_file, output_file, gencode):

    gencode_tb = pysam.TabixFile(gencode)

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Gencode_CDS_start", "Gencode_CDS_end", "Gencode_exon_starts", "Gencode_exon_ends", "Gencode_exon_frames", "Gene_type",
            "Juncmut_hijacked_exon_num","Juncmut_hijacked_SS", "Juncmut_primary_SS", "Juncmut_matching_SS", "SS_shift1", "SS_shift2"
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:
            (juncmut_primary_sj_chr, juncmut_primary_sj_pos) = csvobj["SJ_key"].split(":")
            (juncmut_primary_sj_start, juncmut_primary_sj_end) = list(map(int, juncmut_primary_sj_pos.split('-')))
            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]

            if splice_type in ["Donor+", "Acceptor-"]:
                juncmut_matching_ss = juncmut_primary_sj_end
                juncmut_primary_ss = juncmut_primary_sj_start
            else:
                juncmut_matching_ss = juncmut_primary_sj_start
                juncmut_primary_ss = juncmut_primary_sj_end

            tx_info = define_transcript(splice_type, juncmut_primary_sj_chr, juncmut_matching_ss, gencode_tb, csvobj["Transcript"])
            if tx_info is None: continue
            
            # The ss_shift1 from the corresponding normal splice site. Minus indicates 5' side on the genome. Plus indicates 3'SS.
            ss_shift1 = juncmut_primary_ss - tx_info["Juncmut_hijacked_SS"]
            ss_shift2 = ss_shift1
            # The ss_shift1 from the corresponding normal splice site. Minus indicates the exon side on the pre-mRNA. Minus indicates exon side
            if splice_type in ["Donor-", "Acceptor+"]:
                ss_shift2 = -ss_shift1

            # Is_GTAG_creation
            Is_GTAG_creation = "True"
            if csvobj["Is_GTAG_creation"].startswith("non-"):
                Is_GTAG_creation = "False"

            # gene_type
            gene_type = "Coding"
            if tx_info["Gencode_CDS_start"] == tx_info["Gencode_CDS_end"]:
                gene_type = "Non-coding"

            csvobj["Gencode_CDS_start"] = tx_info["Gencode_CDS_start"]
            csvobj["Gencode_CDS_end"] = tx_info["Gencode_CDS_end"]
            csvobj["Gencode_exon_starts"] = tx_info["Gencode_exon_starts"]
            csvobj["Gencode_exon_ends"] = tx_info["Gencode_exon_ends"]
            csvobj["Gencode_exon_frames"] = tx_info["Gencode_exon_frames"]
            csvobj["Gene_type"] = gene_type
            csvobj["Juncmut_hijacked_exon_num"] = tx_info["Juncmut_hijacked_exon_num"]
            csvobj["Juncmut_hijacked_SS"] = tx_info["Juncmut_hijacked_SS"]
            csvobj["Juncmut_primary_SS"] = juncmut_primary_ss
            csvobj["Juncmut_matching_SS"] = juncmut_matching_ss
            csvobj["SS_shift1"] = ss_shift1
            csvobj["SS_shift2"] = ss_shift2
            csvwriter.writerow(csvobj)

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gencode = sys.argv[3]
    #mane = sys.argv[4]

    sjclass_transcript(input_file, output_file, gencode)
