#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_assadj(input_file, output_file):
    import csv
    print(input_file)
    with open(input_file) as hin, open(output_file, 'w') as hout:
        # [TODO] set sample
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=[
            "SJ_key", "Start_ori", "End_ori", "Created_motif", "SJ_strand", "SJ_read_count", 
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
            
            if "e" in csvobj["Is_Boundary_1"] or "s" in csvobj["Is_Boundary_1"]:
                offset = csvobj["Offset_1"].split(';', 2)[0]
            elif "s" in csvobj["Is_Boundary_2"] or "e" in csvobj["Is_Boundary_2"]:
                offset = csvobj["Offset_2"].split(';',2)[0]
            else:
                raise Exception("juncmut_assadj.py: Unexpected data format")

            if offset == "*":
                print("juncmut_assadj.py: Unexpected data format, offset == '*'")
                offset = 0
            
            if "e" in csvobj["Is_Boundary_1"] or "s" in csvobj["Is_Boundary_2"]:
                strand = "+"
            elif "s" in csvobj["Is_Boundary_1"] or "e" in csvobj["Is_Boundary_2"]:
                strand = "-"
            else:
                continue

            out_csvobj = {
                "SJ_key": "%s,%d,%d" % (csvobj["SJ_1"], int(csvobj["SJ_2"]) - 1 * int(offset), int(csvobj["SJ_3"]) - 1 * int(offset)),
                "Start_ori": csvobj["SJ_2"],
                "End_ori": csvobj["SJ_3"],
                "Created_motif": created_motif,
                "SJ_read_count": csvobj["SJ_7"],
                "SJ_strand": strand
            }
            csvwriter.writerow(out_csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    juncmut_assadj(input_file, output_file)
