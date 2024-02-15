#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_freq(input_file, output_file, original_sj_file, read_num_thres, freq_thres):
    import csv
    
    junc_obj = {}
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            sj_pos = "%s,%s,%s" % (csvobj["Chr"], csvobj["SJ_2"], csvobj["SJ_3"])
            splicing_type = "%s,%s" % (csvobj["Splicing_class"], csvobj["SJ_strand"])
            if not sj_pos in junc_obj:
                junc_obj[sj_pos] = {
                    "SJ_read_count": 0
                }
            if not splicing_type in junc_obj[sj_pos]:
                junc_obj[sj_pos][splicing_type] = {
                    "Start_ori": [],
                    "End_ori": [],
                    "Sample": [],
                }
            junc_obj[sj_pos][splicing_type]["Start_ori"].extend(csvobj["Start_ori"].split(";"))
            junc_obj[sj_pos][splicing_type]["End_ori"].extend(csvobj["End_ori"].split(";"))
            junc_obj[sj_pos][splicing_type]["Sample"].append(csvobj["Sample"])
            junc_obj[sj_pos]["SJ_read_count"] += int(csvobj["SJ_read_count"])
    
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
            "Chr", "Start", "End", "Start_list", "End_list", "Sample", "Splicing_class", "SJ_strand", "SJ_read_count", "SJ_depth", "SJ_freq", 
        ])
        csvwriter.writeheader()

        for sj_pos in sorted(list(junc_obj.keys())):
            # [TODO] fix this case
            if len(junc_obj[sj_pos]) > 2:
                raise Exception("juncmut_freq.py: Unexpected data %s key:%s" % (input_file, sj_pos))

            (chr, sj_2, sj_3) = sj_pos.split(",")
            sj_read_count = junc_obj[sj_pos]["SJ_read_count"]

            for key in sorted(list(junc_obj[sj_pos].keys())):
                if key == 'SJ_read_count': continue
                splicing_type = key
                (splicing_class, sj_strand) = splicing_type.split(",")

                start_ls = sorted(list(set(junc_obj[sj_pos][splicing_type]["Start_ori"]+[sj_2])))
                end_ls = sorted(list(set(junc_obj[sj_pos][splicing_type]["End_ori"]+[sj_3])))

                out_csvobj = {
                    "Chr": chr,
                    "Start": sj_2,
                    "End": sj_3,
                    "Start_list": ";".join(start_ls),
                    "End_list": ";".join(end_ls),
                    "Sample": "",
                    "Splicing_class": splicing_class,
                    "SJ_strand": sj_strand,
                    "SJ_read_count": sj_read_count,
                }

                total = 0
                # strand=+ 5'SS end-side or strand=- 3'SS end-side
                if ("5" in splicing_class and "+" in sj_strand) \
                or ("3" in splicing_class and "-" in sj_strand):
                    for pos in end_ls:
                        total += original_sj_ends.get("%s_%s" % (chr, pos), 0)
                # strand=- 5'SS start-side or strand=+ 3'SS start-side
                elif ("5" in splicing_class and "-" in sj_strand) \
                or ("3" in splicing_class and "+" in sj_strand):
                    for pos in start_ls:
                        total += original_sj_starts.get("%s_%s" % (chr, pos), 0)
                else:
                    continue

                freq = sj_read_count/total
                if sj_read_count >= read_num_thres and freq >= freq_thres:
                    out_csvobj["SJ_depth"] = total
                    out_csvobj["SJ_freq"] = freq
                    csvwriter.writerow(out_csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    original_sj_file = sys.argv[3]
    juncmut_freq(input_file, output_file, original_sj_file, 3, 0.05)
