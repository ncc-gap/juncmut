#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_freq(input_file, output_file, original_sj_file, read_num_thres, freq_thres):
    import csv
    import gzip

    junc_obj = {}
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            sj_key = csvobj["SJ_key"]
            if not sj_key in junc_obj:
                junc_obj[sj_key] = {
                    "SJ_read_count_total": 0,
                    "Start_ori": [],
                    "End_ori": [],
                    "Created_motif": csvobj["Created_motif"],
                    "SJ_strand": csvobj["SJ_strand"],
                    "Transcript": csvobj["Transcript"],
                    "Gene": csvobj["Gene"],
                    "MANE": csvobj["MANE"],
                }
            junc_obj[sj_key]["Start_ori"].extend(csvobj["Start_ori"].split(";"))
            junc_obj[sj_key]["End_ori"].extend(csvobj["End_ori"].split(";"))
            junc_obj[sj_key]["SJ_read_count_total"] += int(csvobj["SJ_read_count"])
    
    original_sj_starts = {}
    original_sj_ends = {}
    with open(original_sj_file) as hin_sj:
        for row in hin_sj:
            F = row.rstrip("\n").split("\t")
            chrom = F[0]
            start = F[1]
            end = F[2]
            read = int(F[6])
            start_pos = "%s_%s" % (chrom, start)
            if not start_pos in original_sj_starts:
                original_sj_starts[start_pos] = read
            else:
                original_sj_starts[start_pos] += read
            
            end_pos = "%s_%s" % (chrom, end)
            if not end_pos in original_sj_ends:
                original_sj_ends[end_pos] = read
            else:
                original_sj_ends[end_pos] += read

    with open(output_file, 'w') as hout:
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=[
            "SJ_key", "Created_motif", "SJ_strand", "Transcript", "Gene", "MANE", "SJ_read_count", "SJ_depth", "SJ_freq"
        ])
        csvwriter.writeheader()

        for sj_key in sorted(list(junc_obj.keys())):
            (sj_key_chr, sj_key_start, sj_key_end) = sj_key.split(",")

            depth = 0
            splice_type = junc_obj[sj_key]["Created_motif"] + junc_obj[sj_key]["SJ_strand"]
            # strand=+ 5'SS end-side or strand=- 3'SS end-side
            if splice_type == "Donor+" or splice_type == "Acceptor-":
                end_list = sorted(list(set(junc_obj[sj_key]["End_ori"]+[sj_key_end])))
                for pos in end_list:
                    depth += original_sj_ends.get("%s_%s" % (sj_key_chr, pos), 0)
            # strand=- 5'SS start-side or strand=+ 3'SS start-side
            else:
                start_list = sorted(list(set(junc_obj[sj_key]["Start_ori"]+[sj_key_start])))
                for pos in start_list:
                    depth += original_sj_starts.get("%s_%s" % (sj_key_chr, pos), 0)

            freq = junc_obj[sj_key]["SJ_read_count_total"]/depth

            if junc_obj[sj_key]["SJ_read_count_total"] >= read_num_thres and freq >= freq_thres:
                out_csvobj = {
                    "SJ_key": sj_key,
                    "Created_motif": junc_obj[sj_key]["Created_motif"],
                    "SJ_strand": junc_obj[sj_key]["SJ_strand"],
                    "SJ_read_count": junc_obj[sj_key]["SJ_read_count_total"],
                    "SJ_depth": depth,
                    "SJ_freq": freq,
                    "Transcript": junc_obj[sj_key]["Transcript"],
                    "Gene": junc_obj[sj_key]["Gene"],
                    "MANE": junc_obj[sj_key]["MANE"],
                }
                csvwriter.writerow(out_csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    original_sj_file = sys.argv[3]

    juncmut_freq(input_file, output_file, original_sj_file, 3, 0.05)
