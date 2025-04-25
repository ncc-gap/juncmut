#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_assadj(input_file, output_file):
    import csv
    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=[
            "SJ_key", "Start_ori", "End_ori", "Created_motif", "SJ_strand", "Transcript", "Gene", "MANE", "SJ_read_count"
        ])
        csvwriter.writeheader()
        
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            if not "alternative" in csvobj["Splicing_Class"].lower():
                continue

            created_motif = ""
            if "5'SS" in csvobj["Splicing_Class"]:
                created_motif = "Donor"
            elif "3'SS" in csvobj["Splicing_Class"]:
                created_motif = "Acceptor"
            else:
                continue
            
            # intron pos1
            index_1 = None
            mane_select_list_1 = csvobj["MANE_1"].split(';')
            is_boundary_list_1 = csvobj["Is_Boundary_1"].split(';')

            if "MANE_Select" in mane_select_list_1:
                for i,is_boundary in enumerate(is_boundary_list_1):
                    if mane_select_list_1[i] == "MANE_Select" and is_boundary in ["s", "e"]:
                        index_1 = i
                        mane = "MANE_Select"
                        break
            else:
                for i,is_boundary in enumerate(is_boundary_list_1):
                    if is_boundary in ["s", "e"]:
                        index_1 = i
                        mane = "NA"
                        break

            if not index_1 is None:
                transcript = csvobj["Gene_1"].split(';')[index_1]
                symbol = csvobj["Symbol_1"].split(';')[index_1]
                offset = csvobj["Offset_1"].split(';')[index_1]
                if is_boundary_list_1[index_1] == "s":
                    strand = "-"
                else:
                    strand = "+"

            # intron pos2
            index_2 = None
            if index_1 is None:
                mane_select_list_2 = csvobj["MANE_2"].split(';')
                is_boundary_list_2 = csvobj["Is_Boundary_2"].split(';')

                if "MANE_Select" in mane_select_list_2:
                    for i,is_boundary in enumerate(is_boundary_list_2):
                        if mane_select_list_2[i] == "MANE_Select" and is_boundary in ["s", "e"]:
                            index_2 = i
                            mane = "MANE_Select"
                            break
                else:
                    for i,is_boundary in enumerate(is_boundary_list_2):
                        if is_boundary in ["s", "e"]:
                            index_2 = i
                            mane = "NA"
                            break

                if not index_2 is None:
                    transcript = csvobj["Gene_2"].split(';')[index_2]
                    symbol = csvobj["Symbol_2"].split(';')[index_2]
                    offset = csvobj["Offset_2"].split(';')[index_2]
                    if csvobj["Is_Boundary_2"].split(';')[index_2] == "s":
                        strand = "+"
                    else:
                        strand = "-"

            if index_1 is None and index_2 is None:
                print(csvobj)
                raise Exception("juncmut_assadj.py: Unexpected data format")

            out_csvobj = {
                "SJ_key": "%s,%d,%d" % (csvobj["SJ_1"], int(csvobj["SJ_2"]) - 1 * int(offset), int(csvobj["SJ_3"]) - 1 * int(offset)),
                "Start_ori": csvobj["SJ_2"],
                "End_ori": csvobj["SJ_3"],
                "Created_motif": created_motif,
                "SJ_read_count": csvobj["SJ_7"],
                "SJ_strand": strand,
                "Transcript": transcript,
                "Gene": symbol,
                "MANE": mane,
            }
            csvwriter.writerow(out_csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    juncmut_assadj(input_file, output_file)
