#!/usr/bin/env python3

import csv

def annot_cgc(input_file, output_file, cgc_file):
    
    #make dictionary
    cgc2role = {}
    with open(cgc_file, 'r') as hin:
        dreader = csv.DictReader(hin, delimiter = '\t')
        for F in dreader:
 
            if F["Tier"] != '1': continue

            mutation_flag = False
            if 'F' in F["Mutation Types"]: mutation_flag = True
            if 'Mis' in F["Mutation Types"]: mutation_flag = True
            if 'N' in F["Mutation Types"]: mutation_flag = True 
            if 'S' in F["Mutation Types"]: mutation_flag = True  
            if not mutation_flag: continue

            role_value = "-"
            if "oncogene" in F["Role in Cancer"] and "TSG" in F["Role in Cancer"]:
                role_value = "oncogene,TSG"
            elif "oncogene" in F["Role in Cancer"]:
                role_value = "oncogene"
            elif "TSG" in F["Role in Cancer"]:
                role_value = "TSG"

            cgc2role[F["Gene Symbol"]] = role_value

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', quotechar='"', fieldnames=csvreader.fieldnames + ["Is_CGC", "CGC_role"])
        csvwriter.writeheader()

        for csvobj in csvreader:
            # [TODO] Requires checking
            #csvobj["Is_CGC"] = "TRUE" if csvobj["Gencode_name2"] in cgc2role else "FALSE"
            #csvobj["CGC_role"] = cgc2role[csvobj["Gencode_name2"]] if csvobj["Gencode_name2"] in cgc2role else "-"
            csvobj["Is_CGC"] = "TRUE" if csvobj["Gene"] in cgc2role else "FALSE"
            csvobj["CGC_role"] = cgc2role[csvobj["Gene"]] if csvobj["Gene"] in cgc2role else "---"
            csvwriter.writerow(csvobj)
 
if __name__== "__main__":
    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    cgc_file = sys.argv[3]
    
    annot_cgc(input_file, output_file, cgc_file)
