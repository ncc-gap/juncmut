#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_freq(input_file, output_file, original_sj_file, read_num_thres, freq_thres):
    import csv
    
    junc_obj = {}
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            sj_key = csvobj["SJ_key"]
            splicing_type = csvobj["Created_motif"] + csvobj["SJ_strand"]
            if not sj_key in junc_obj:
                junc_obj[sj_key] = {
                    "SJ_read_count_total": 0
                }
            if not splicing_type in junc_obj[sj_key]:
                junc_obj[sj_key][splicing_type] = {
                    "Start_ori": [],
                    "End_ori": [],
                }
            junc_obj[sj_key][splicing_type]["Start_ori"].extend(csvobj["Start_ori"].split(";"))
            junc_obj[sj_key][splicing_type]["End_ori"].extend(csvobj["End_ori"].split(";"))
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
            "SJ_key", "Created_motif", "SJ_strand", "SJ_read_count", "SJ_depth", "SJ_freq", 
        ])
        csvwriter.writeheader()

        for sj_key in sorted(list(junc_obj.keys())):
            # [TODO] fix this case
            if len(junc_obj[sj_key]) > 2:
                raise Exception("juncmut_freq.py: Unexpected data %s key:%s" % (input_file, sj_key))

            (sj_key_chr, sj_key_start, sj_key_end) = sj_key.split(",")
            sj_read_count = junc_obj[sj_key]["SJ_read_count_total"]

            for key in sorted(list(junc_obj[sj_key].keys())):
                if key == 'SJ_read_count_total':
                    continue
                splicing_type = key
                start_ls = sorted(list(set(junc_obj[sj_key][splicing_type]["Start_ori"]+[sj_key_start])))
                end_ls = sorted(list(set(junc_obj[sj_key][splicing_type]["End_ori"]+[sj_key_end])))

                depth = 0
                # strand=+ 5'SS end-side or strand=- 3'SS end-side
                if splicing_type == "Donor+" or splicing_type == "Acceptor-":
                    for pos in end_ls:
                        depth += original_sj_ends.get("%s_%s" % (sj_key_chr, pos), 0)
                # strand=- 5'SS start-side or strand=+ 3'SS start-side
                else:
                    for pos in start_ls:
                        depth += original_sj_starts.get("%s_%s" % (sj_key_chr, pos), 0)

                freq = sj_read_count/depth
                if sj_read_count >= read_num_thres and freq >= freq_thres:
                    out_csvobj = {
                        "SJ_key": sj_key,
                        "Created_motif": splicing_type[0:-1],
                        "SJ_strand": splicing_type[-1],
                        "SJ_read_count": sj_read_count,
                        "SJ_depth": depth,
                        "SJ_freq": freq,
                    }
                    csvwriter.writerow(out_csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    original_sj_file = sys.argv[3]
    juncmut_freq(input_file, output_file, original_sj_file, 3, 0.05)
