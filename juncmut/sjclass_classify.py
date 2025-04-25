#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 09:01:08 2021

@author: Naoko Iida
"""
import subprocess
import pysam
import csv
import os

SJ_HEADER = [
  'chromosome',
  'first_base',
  'last_base',
  'strand',
  'intron_motif',
  'annotated',
  'number_of_uniquely_mapping_reads',
  'number_of_multi-mapping_reads',
  'maximum_spliced_alignment',
]

def sjclass_classify(input_file, output_file, bam_file, sj_file, depth_th):

    with open(output_file + ".SJ.out.bed", 'w') as hout_bed:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", sj_file], stdout = hout_bed)
    with open(output_file + ".SJ.out.bed.gz", 'wb') as hout_bedgz:
        subprocess.check_call(["bgzip", "-f", "-c", output_file + ".SJ.out.bed"], stdout = hout_bedgz)

    subprocess.check_call(["tabix", "-p", "bed", output_file + ".SJ.out.bed.gz"])
    subprocess.check_call(["rm", "-rf", output_file + ".SJ.out.bed"])

    sj_tb = pysam.TabixFile(output_file + ".SJ.out.bed.gz")
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Juncmut_primary_SJ", "Juncmut_hijacked_SJ", "Juncmut_secondary_SS", "Juncmut_secondary_SJ", 
            "Closed_exon_num", "Juncmut_predicted_splicing_type", "Num_skipped_exon",
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:

            juncmut_predicted_splicing_type = "NA"
            juncmut_secondary_ss = "NA"
            juncmut_secondary_sj = "NA"
            num_skipped_exon = "NA"
            closed_exon_num = "NA"

            csvobj["Juncmut_primary_SJ"] = "NA"
            csvobj["Juncmut_hijacked_SJ"] = "NA"
            csvobj["Juncmut_secondary_SS"] = juncmut_secondary_ss
            csvobj["Juncmut_secondary_SJ"] = juncmut_secondary_sj
            csvobj["Closed_exon_num"] = closed_exon_num
            csvobj["Juncmut_predicted_splicing_type"] = juncmut_predicted_splicing_type
            csvobj["Num_skipped_exon"] = num_skipped_exon
            if csvobj["Gencode_exon_starts"] == "NA":
                csvwriter.writerow(csvobj)
                continue

            (juncmut_primary_sj_chr, juncmut_primary_sj_pos) = csvobj["SJ_key"].split(":")
            (juncmut_primary_sj_start, juncmut_primary_sj_end) = list(map(int, juncmut_primary_sj_pos.split('-')))
            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]
            Juncmut_hijacked_ss = int(csvobj["Juncmut_hijacked_SS"])
            juncmut_primary_ss = int(csvobj["Juncmut_primary_SS"])
            juncmut_matching_ss = int(csvobj["Juncmut_matching_SS"])

            gencode_exon_starts = list(map(int, csvobj["Gencode_exon_starts"].rstrip(',').split(',')))
            gencode_exon_ends = list(map(int, csvobj["Gencode_exon_ends"].rstrip(',').split(',')))
            (hijacked_exon_num, gencode_exon_count) = list(map(int, csvobj["Juncmut_hijacked_exon_num"].replace(',','').split('/')))

            # o--->
            if splice_type in ["Donor+", "Acceptor-"]:
                if juncmut_primary_ss - 1 in gencode_exon_ends:
                    csvwriter.writerow(csvobj)
                    continue

                for e_count in range(hijacked_exon_num, -1, -1):
                    if gencode_exon_starts[e_count] + 1 <= juncmut_primary_ss:
                        closed_exon_num = e_count
                        num_skipped_exon = hijacked_exon_num - closed_exon_num
                        break

                if closed_exon_num == "NA":
                    juncmut_predicted_splicing_type = "Ambiguous_outgene"

                elif gencode_exon_starts[closed_exon_num] < juncmut_primary_ss <= gencode_exon_ends[closed_exon_num]:
                    juncmut_predicted_splicing_type = "Shortening_exon"

                else:
                    # cryptic exon
                    sj_records = sj_tb.fetch(region = "%s:%d-%d" % (juncmut_primary_sj_chr, gencode_exon_ends[closed_exon_num], juncmut_primary_ss))
                    if sj_records == None:
                        sj_records = []
                    for record in sj_records:
                        R = record.split('\t')
                        sj_start = int(R[SJ_HEADER.index("first_base")])
                        sj_end = int(R[SJ_HEADER.index("last_base")])
                        sj_read_n = int(R[SJ_HEADER.index("number_of_uniquely_mapping_reads")])
                        dist2sj = juncmut_primary_ss - sj_end
                        # sj_SS is within norm_SS and alt_SS,considering the shift of alignment(2 bp). The size of exon is likely <=300 bp.
                        # exon_e is the last base of the exon
                        if gencode_exon_ends[closed_exon_num] < sj_end < juncmut_primary_ss and dist2sj <= 300 and sj_read_n > 0 \
                        and gencode_exon_ends[closed_exon_num] - 2 <= sj_start <= gencode_exon_ends[closed_exon_num] + 2:
                            juncmut_predicted_splicing_type = "Cryptic_exon"
                            mis_dis = sj_start - (gencode_exon_ends[closed_exon_num] + 1)
                            juncmut_secondary_ss = sj_end + mis_dis
                            juncmut_secondary_sj = "%s:%d-%d" % (juncmut_primary_sj_chr, gencode_exon_ends[closed_exon_num] + 1, juncmut_secondary_ss)
                            break

                if juncmut_predicted_splicing_type == "NA":
                    # depth
                    region = "%s:%d-%d" % (juncmut_primary_sj_chr, gencode_exon_ends[closed_exon_num] + 1, juncmut_primary_ss - 1)
                    with open(output_file + ".depth.txt", 'w') as hout_depth:
                        subprocess.run(["samtools", "depth", "-a", "-r", region, bam_file], stdout = hout_depth)

                    #df = pd.read_csv(output_file + ".depth.txt", sep="\t", header=None)
                    #df_a = df[df.iloc[:,2] >= depth_th]
                    #a = set(df_a.iloc[:,1].tolist())
                    a = []
                    with open(output_file + ".depth.txt") as hin_depth:
                        for row in hin_depth:
                            F = row.rstrip("\n").split("\t")
                            if len(F) < 3:
                                continue
                            depth_pos = int(F[1])
                            depth_value = int(F[2])
                            if depth_value >= depth_th:
                                a.append(depth_pos)
                    c = range(gencode_exon_ends[closed_exon_num] + 1, juncmut_primary_ss, 1)
                    # lengthen exon
                    if set(c) <= set(a):
                        juncmut_predicted_splicing_type = "Lengthening_exon"
                    else:
                        juncmut_predicted_splicing_type = "Ambiguous_termination"

                    os.remove(output_file +".depth.txt")
            
            # <---o
            elif splice_type in ["Donor-", "Acceptor+"]:
                if juncmut_primary_ss in gencode_exon_starts:
                    csvwriter.writerow(csvobj)
                    continue

                for e_count in range(hijacked_exon_num, gencode_exon_count):
                    if gencode_exon_ends[e_count] > juncmut_primary_ss:
                        closed_exon_num = e_count
                        num_skipped_exon = closed_exon_num - hijacked_exon_num
                        break

                if closed_exon_num == "NA":
                    juncmut_predicted_splicing_type = "Ambiguous_outgene"
                    #print("Ambiguous_outgene")
                elif gencode_exon_starts[closed_exon_num] < juncmut_primary_ss < gencode_exon_ends[closed_exon_num]:
                    juncmut_predicted_splicing_type = "Shortening_exon"
                    #print("Shortening_exon")
                else:
                    #cryptic exon
                    sj_records = sj_tb.fetch(region = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_primary_ss, gencode_exon_starts[closed_exon_num]))
                    if sj_records == None:
                        sj_records == []
                    for record in sj_records:
                        R = record.split('\t')
                        sj_start = int(R[SJ_HEADER.index("first_base")])
                        sj_end = int(R[SJ_HEADER.index("last_base")])
                        sj_read_n = int(R[SJ_HEADER.index("number_of_uniquely_mapping_reads")])
                        dist2sj = sj_start - juncmut_primary_ss
                        # sj_SS is within norm_SS and alt_SS,considering the shift of alignment(2 bp). The size of exon is likely <=300 bp.
                        if juncmut_primary_ss+1 < sj_start < gencode_exon_starts[closed_exon_num] and dist2sj <= 300 and sj_read_n > 0 \
                        and gencode_exon_starts[closed_exon_num] - 2 <= sj_end <= gencode_exon_starts[closed_exon_num] + 2:
                            juncmut_predicted_splicing_type = "Cryptic_exon"
                            mis_dis = sj_end - gencode_exon_starts[closed_exon_num]
                            juncmut_secondary_ss = sj_start + mis_dis
                            juncmut_secondary_sj = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_secondary_ss, gencode_exon_starts[closed_exon_num])
                            break
                # depth
                if juncmut_predicted_splicing_type == "NA":
                    region = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_primary_ss + 1, gencode_exon_starts[closed_exon_num])
                    with open(output_file + ".depth.txt", 'w') as hout_depth:
                        subprocess.run(["samtools", "depth", "-a", "-r", region, bam_file], stdout = hout_depth)

                    #df = pd.read_csv(output_file + ".depth.txt", sep="\t", header=None)
                    #df_a = df[df.iloc[:,2] >= depth_th]
                    #a = set(df_a.iloc[:,1].tolist())
                    a = []
                    with open(output_file + ".depth.txt") as hin_depth:
                        for row in hin_depth:
                            F = row.rstrip("\n").split("\t")
                            if len(F) < 3:
                                continue
                            depth_pos = int(F[1])
                            depth_value = int(F[2])
                            if depth_value >= depth_th:
                                a.append(depth_pos)
                    c = range(juncmut_primary_ss + 1, gencode_exon_starts[closed_exon_num] + 1, 1)
                    # lengthen exon
                    if set(c) <= set(a):
                        juncmut_predicted_splicing_type = "Lengthening_exon"
                    else:
                        juncmut_predicted_splicing_type = "Ambiguous_termination"

                    os.remove(output_file +".depth.txt")
            


            # SJ
            if splice_type in ["Donor+", "Acceptor-"]:
                csvobj["Juncmut_primary_SJ"] = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_primary_ss, juncmut_matching_ss)
                csvobj["Juncmut_hijacked_SJ"] = "%s:%d-%d" % (juncmut_primary_sj_chr, Juncmut_hijacked_ss, juncmut_matching_ss)

            elif splice_type in ["Donor-", "Acceptor+"]:
                csvobj["Juncmut_primary_SJ"] = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_matching_ss, juncmut_primary_ss)
                csvobj["Juncmut_hijacked_SJ"] = "%s:%d-%d" % (juncmut_primary_sj_chr, juncmut_matching_ss, Juncmut_hijacked_ss)

            csvobj["Juncmut_secondary_SS"] = juncmut_secondary_ss
            csvobj["Juncmut_secondary_SJ"] = juncmut_secondary_sj

            # Closed_exon_num
            if closed_exon_num == "NA":
                csvobj["Closed_exon_num"] = closed_exon_num
            else:
                csvobj["Closed_exon_num"] = "%d/%d," % (closed_exon_num, gencode_exon_count)

            csvobj["Juncmut_predicted_splicing_type"] = juncmut_predicted_splicing_type
            csvobj["Num_skipped_exon"] = num_skipped_exon

            csvwriter.writerow(csvobj)

    sj_tb.close()
    subprocess.check_call(["rm", "-rf", output_file + ".SJ.out.bed.gz"])
    subprocess.check_call(["rm", "-rf", output_file + ".SJ.out.bed.gz.tbi"])

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bam_file = sys.argv[3]
    sj_file = sys.argv[4]

    sjclass_classify(input_file, output_file, bam_file, sj_file, 1)
