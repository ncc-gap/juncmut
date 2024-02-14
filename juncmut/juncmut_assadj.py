#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_assadj(input_file, output_file):
    import csv
    
    with open(input_file) as hin, open(output_file, 'w') as hout:
        # [TODO] set sample
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=[
            "Chr", "Start_ori", "End_ori", "SJ_2", "SJ_3", "Sample", "Splicing_class", "SJ_strand", "SJ_read_count", 
        ])
        csvwriter.writeheader()
        
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            if not "alternative" in csvobj["Splicing_Class"].lower():
                continue
            
            out_csvobj = {
                "Chr": csvobj["SJ_1"],
                "SJ_2": csvobj["SJ_2"],
                "SJ_3": csvobj["SJ_3"],
                "Sample": "",
                "Splicing_class": csvobj["Splicing_Class"],
                "SJ_read_count": csvobj["SJ_7"],
            }
            
            if "e" in csvobj["Is_Boundary_1"] or "s" in csvobj["Is_Boundary_1"]:
                offset = csvobj["Offset_1"].split(';', 2)[0]
            elif "s" in csvobj["Is_Boundary_2"] or "e" in csvobj["Is_Boundary_2"]:
                offset = csvobj["Offset_2"].split(';',2)[0]
            else:
                raise Exception("juncmut_assadj.py: Unexpected data format")
            
            if "e" in csvobj["Is_Boundary_1"] or "s" in csvobj["Is_Boundary_2"]:
                strand = "+"
            elif "s" in csvobj["Is_Boundary_1"] or "e" in csvobj["Is_Boundary_2"]:
                strand = "-"
            
            out_csvobj["Start_ori"] = int(csvobj["SJ_2"]) - 1 * int(offset)
            out_csvobj["End_ori"] = int(csvobj["SJ_3"]) - 1 * int(offset)
            out_csvobj["SJ_strand"] = strand
            csvwriter.writerow(out_csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    juncmut_assadj(input_file, output_file)
