#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:01:08 2021

@author: Naoko Iida
"""

import json
import gzip
import csv

def min_cds_edge(start_end_codon, cds_list_edge, default_value):

    if start_end_codon is None and cds_list_edge is None:
        return default_value

    if start_end_codon is None:
        return cds_list_edge

    if cds_list_edge is None:
        return start_end_codon

    return str(min(int(start_end_codon), int(cds_list_edge)))

def max_cds_edge(start_end_codon, cds_list_edge, default_value):

    if start_end_codon is None and cds_list_edge is None:
        return default_value

    if start_end_codon is None:
        return cds_list_edge

    if cds_list_edge is None:
        return start_end_codon

    return str(max(int(start_end_codon), int(cds_list_edge)))

def load_transcript(input_file, gencode_gene_file):
    if gencode_gene_file.endswith("gtf.gz"):
        key_value_split = " "
    elif gencode_gene_file.endswith("gff3.gz"):
        key_value_split = "="
    else:
        raise Exception("unexcepted gencode gene file")

    target_transcripts = []
    with open(input_file, 'r') as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            if not csvobj["Transcript"] in target_transcripts:
                target_transcripts.append(csvobj["Transcript"])

    transcripts = {}
    with gzip.open(gencode_gene_file, 'rt') as hin:

        prev_transcript_id = ""
        start_codon_start = None
        start_codon_end = None
        stop_codon_start = None
        stop_codon_end = None
        cds_starts = []
        cds_ends = []
        exon_starts = []
        exon_ends = []
        strand = ""
        for line in hin:
            if line.startswith("#"):
                continue

            F = line.rstrip('\n').split('\t')
            feature_type = F[2]
            if not feature_type in ["exon", "start_codon", "stop_codon", "CDS"]:
                continue

            transcript_id = ""
            for item in F[8].rstrip(";").split(";"):
                (key, value) = item.strip(" ").replace('"', '').split(key_value_split)
                if key == "transcript_id":
                    transcript_id = value
                    break

            if prev_transcript_id != "" and prev_transcript_id != transcript_id:
                if strand == "+":
                    transcripts[prev_transcript_id] = {
                        "Gencode_CDS_start": min_cds_edge(start_codon_start, cds_starts[0] if cds_starts else None, exon_starts[0]), 
                        "Gencode_CDS_end": max_cds_edge(stop_codon_end, cds_ends[-1] if cds_ends else None, exon_starts[0]), 
                        "Gencode_exon_count": len(exon_starts), 
                        "Gencode_exon_starts": ",".join(exon_starts) + ",", 
                        "Gencode_exon_ends": ",".join(exon_ends) + ",",
                    }

                elif strand == "-":
                    transcripts[prev_transcript_id] = {
                        "Gencode_CDS_start": min_cds_edge(stop_codon_start, cds_starts[-1] if cds_starts else None, exon_starts[-1]), 
                        "Gencode_CDS_end": max_cds_edge(start_codon_end, cds_ends[0] if cds_ends else None, exon_starts[-1]), 
                        "Gencode_exon_count": len(exon_starts), 
                        "Gencode_exon_starts": ",".join(reversed(exon_starts)) + ",", 
                        "Gencode_exon_ends": ",".join(reversed(exon_ends)) + ",",
                    }

                start_codon_start = None
                start_codon_end = None
                stop_codon_start = None
                stop_codon_end = None
                cds_starts = []
                cds_ends = []
                exon_starts = []
                exon_ends = []
                strand = ""
                prev_transcript_id = ""
                if len(transcripts) == len(target_transcripts):
                    break

            if not transcript_id in target_transcripts:
                continue

            prev_transcript_id = transcript_id
            start = str(int(F[3]) - 1)
            end = F[4]
            strand = F[6]
            if feature_type == "start_codon":
                start_codon_start = start
                start_codon_end = end
            elif feature_type == "stop_codon":
                stop_codon_start = start
                stop_codon_end = end
            elif feature_type == "CDS":
                cds_starts.append(start)
                cds_ends.append(end)
            else:
                exon_starts.append(start)
                exon_ends.append(end)

        if prev_transcript_id != "":
            if strand == "+":
                transcripts[prev_transcript_id] = {
                    "Gencode_CDS_start": min_cds_edge(start_codon_start, cds_starts[0] if cds_starts else None, exon_starts[0]), 
                    "Gencode_CDS_end": max_cds_edge(stop_codon_end, cds_ends[-1] if cds_ends else None, exon_starts[0]), 
                    "Gencode_exon_count": len(exon_starts), 
                    "Gencode_exon_starts": ",".join(exon_starts) + ",", 
                    "Gencode_exon_ends": ",".join(exon_ends) + ",",
                }

            elif strand == "-":
                transcripts[prev_transcript_id] = {
                    "Gencode_CDS_start": min_cds_edge(stop_codon_start, cds_starts[-1] if cds_starts else None, exon_starts[-1]), 
                    "Gencode_CDS_end": max_cds_edge(start_codon_end, cds_ends[0] if cds_ends else None, exon_starts[-1]), 
                    "Gencode_exon_count": len(exon_starts), 
                    "Gencode_exon_starts": ",".join(reversed(exon_starts)) + ",", 
                    "Gencode_exon_ends": ",".join(reversed(exon_ends)) + ",",
                }

    return transcripts

def sjclass_transcript(input_file, output_file, gencode_gene_file):

    transcripts = load_transcript(input_file, gencode_gene_file)

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Gencode_CDS_start", "Gencode_CDS_end", "Gencode_exon_starts", "Gencode_exon_ends", "Gene_type",
            "Juncmut_hijacked_exon_num","Juncmut_hijacked_SS", "Juncmut_primary_SS", "Juncmut_matching_SS", "SS_shift1", "SS_shift2",
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

            if not csvobj["Transcript"] in transcripts:
                csvobj["Gencode_CDS_start"] = "NA"
                csvobj["Gencode_CDS_end"] = "NA"
                csvobj["Gencode_exon_starts"] = "NA"
                csvobj["Gencode_exon_ends"] = "NA"
                csvobj["Gene_type"] = "NA"
                csvobj["Juncmut_hijacked_exon_num"] = "NA"
                csvobj["Juncmut_hijacked_SS"] = "NA"
                csvobj["Juncmut_primary_SS"] = "NA"
                csvobj["Juncmut_matching_SS"] = "NA"
                csvobj["SS_shift1"] = "NA"
                csvobj["SS_shift2"] = "NA"
                csvwriter.writerow(csvobj)
                continue

            gencode_exon_starts = list(map(int, transcripts[csvobj["Transcript"]]["Gencode_exon_starts"].rstrip(',').split(',')))
            gencode_exon_ends = list(map(int, transcripts[csvobj["Transcript"]]["Gencode_exon_ends"].rstrip(',').split(',')))

            juncmut_hijacked_ss = -1
            juncmut_hijacked_exon_num = -1
            if splice_type == "Donor+":
                index_start = gencode_exon_starts.index(juncmut_matching_ss)
                if index_start == 0:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_ends[index_start - 1] + 1
                    # be careful, the index_start is a 0-start count. juncmut_hijacked_exon_num is a 1-start count
                    juncmut_hijacked_exon_num = index_start - 1
            
            elif splice_type == "Donor-":
                index_end = gencode_exon_ends.index(juncmut_matching_ss - 1)
                if index_end == len(gencode_exon_ends) - 1:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_starts[index_end + 1]
                    juncmut_hijacked_exon_num = index_end + 1
            
            elif splice_type == "Acceptor+":
                index_end = gencode_exon_ends.index(juncmut_matching_ss - 1)
                if index_end == len(gencode_exon_ends) - 1:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_starts[index_end + 1]
                    juncmut_hijacked_exon_num = index_end + 1
            
            elif splice_type == "Acceptor-":
                index_start = gencode_exon_starts.index(juncmut_matching_ss)
                if index_start == 0:
                    continue
                else:
                    juncmut_hijacked_ss = gencode_exon_ends[index_start - 1] + 1
                    juncmut_hijacked_exon_num = index_start - 1
            
            if juncmut_hijacked_ss == -1:
                continue

            # The ss_shift1 from the corresponding normal splice site. Minus indicates 5' side on the genome. Plus indicates 3'SS.
            ss_shift1 = juncmut_primary_ss - juncmut_hijacked_ss
            ss_shift2 = ss_shift1
            # The ss_shift1 from the corresponding normal splice site. Minus indicates the exon side on the pre-mRNA. Minus indicates exon side
            if splice_type in ["Donor-", "Acceptor+"]:
                ss_shift2 = -ss_shift1

            # gene_type
            gene_type = "Coding"
            if transcripts[csvobj["Transcript"]]["Gencode_CDS_start"] == transcripts[csvobj["Transcript"]]["Gencode_CDS_end"]:
                gene_type = "Non-coding"

            csvobj["Gencode_CDS_start"] = transcripts[csvobj["Transcript"]]["Gencode_CDS_start"]
            csvobj["Gencode_CDS_end"] = transcripts[csvobj["Transcript"]]["Gencode_CDS_end"]
            csvobj["Gencode_exon_starts"] = transcripts[csvobj["Transcript"]]["Gencode_exon_starts"]
            csvobj["Gencode_exon_ends"] = transcripts[csvobj["Transcript"]]["Gencode_exon_ends"]
            csvobj["Gene_type"] = gene_type
            csvobj["Juncmut_hijacked_exon_num"] = "%d/%s," % (juncmut_hijacked_exon_num, transcripts[csvobj["Transcript"]]["Gencode_exon_count"])
            csvobj["Juncmut_hijacked_SS"] = juncmut_hijacked_ss
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

    sjclass_transcript(input_file, output_file, gencode)
