#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_mutpre(input_file, output_file, reference):
    import csv
    import pysam
    
    ref = pysam.FastaFile(reference) # 0-coordinate
    
    dict5p = {0: 'A', 1: 'G', 2: 'G', 3: 'T', 4: 'AG', 5: 'A', 6: 'G', 7: 'T'}
    dict3p = {0: 'CT', 1: 'CT', 2: 'ATGC', 3: 'CT', 4: 'A', 5: 'G', 6: 'AG'}
    dict5m = {0: 'A', 1: 'C', 2: 'T', 3: 'CT', 4: 'A', 5: 'C', 6: 'C', 7: 'T'}
    dict3m = {0: 'CT', 1: 'C', 2: 'T', 3: 'AG', 4: 'TCGA', 5: 'AG', 6: 'AG'}
    
    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=csvreader.fieldnames + [
            "Ref_motif", "Possive_alt_motif", "Possive_alt_key", "Is_GTAG_creation"
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:
            (sj_key_chr, sj_key_start, sj_key_end) = csvobj["SJ_key"].split(",")
            start = int(sj_key_start)
            end = int(sj_key_end)
            splice_type = csvobj["Created_motif"] + csvobj["SJ_strand"]

            possive_alt_motif = []
            possive_alt_key = []
            is_gtag_creation = False

            if splice_type == "Donor+": #o-->
                ref_motif = ref.fetch(sj_key_chr, start-1-2, start+5).upper()

                temp_possive_alt_key = {}
                for i in range(0,len(ref_motif)):
                    if i == 4 and ref_motif[i] == "G":
                        possive_alt_motif.append("A")
                        temp_possive_alt_key[i] = "%d:G:A" % (start+i-2)

                    elif ref_motif[i] in dict5p[i]:
                        possive_alt_motif.append(".")
                    else:
                        if dict5p[i] == "AG":
                            possive_alt_motif.append("R")
                        else:
                            possive_alt_motif.append(dict5p[i])

                        temp_possive_alt_key[i] = "%d:%s:%s" % (start+i-2, ref_motif[i], dict5p[i])

                if possive_alt_motif[2] == "G" or possive_alt_motif[3] == "T":
                    is_gtag_creation = True
                    if possive_alt_motif[2] == "G":
                        possive_alt_key.append(temp_possive_alt_key[2])

                    if possive_alt_motif[3] == "T":
                        possive_alt_key.append(temp_possive_alt_key[3])

                elif list(set(possive_alt_motif)) == ["."]:
                    possive_alt_key = ['-']
                else:
                    for i in sorted(temp_possive_alt_key.keys()):
                        possive_alt_key.append(temp_possive_alt_key[i])

            elif splice_type == "Donor-": #<--o
                ref_motif = ref.fetch(sj_key_chr, end-1-5, end+2).upper()

                temp_possive_alt_key = {}
                for i in range(0,len(ref_motif)):
                    if ref_motif[i] in dict5m[i]:
                        possive_alt_motif.append(".")
                    else:
                        if dict5m[i] == "CT":
                            possive_alt_motif.append("Y")
                        else:
                            possive_alt_motif.append(dict5m[i])
                        
                        temp_possive_alt_key[i] = "%d:%s:%s" % (end+i-5, ref_motif[i], dict5m[i])

                if possive_alt_motif[4] == "A" or possive_alt_motif[5] == "C":
                    is_gtag_creation = True
                    if possive_alt_motif[4] == "A":
                        possive_alt_key.append(temp_possive_alt_key[4])

                    if possive_alt_motif[5] == "C":
                        possive_alt_key.append(temp_possive_alt_key[5])

                elif list(set(possive_alt_motif)) == ["."]:
                    possive_alt_key = ['-']
                else:
                    for i in sorted(temp_possive_alt_key.keys()):
                        possive_alt_key.append(temp_possive_alt_key[i])
            
            elif splice_type == "Acceptor+": #-->o
                ref_motif = ref.fetch(sj_key_chr, end-1-5, end+1).upper()

                temp_possive_alt_key = {}
                for i in range(0,len(ref_motif)):
                    if ref_motif[i] in dict3p[i]:
                        possive_alt_motif.append(".")
                    else:
                        if dict3p[i] == "AG":
                            possive_alt_motif.append("R")
                        elif dict3p[i] == "CT":
                            possive_alt_motif.append("Y")
                        else:
                            possive_alt_motif.append(dict3p[i])

                        temp_possive_alt_key[i] = "%d:%s:%s" % (end+i-5, ref_motif[i], dict3p[i])

                if possive_alt_motif[4] == "A" or possive_alt_motif[5] == "G":
                    is_gtag_creation = True
                    if possive_alt_motif[4] == "A":
                        possive_alt_key.append(temp_possive_alt_key[4])

                    if possive_alt_motif[5] == "G":
                        possive_alt_key.append(temp_possive_alt_key[5])

                elif list(set(possive_alt_motif)) == ["."]:
                    possive_alt_key = ['-']
                else:
                    for i in sorted(temp_possive_alt_key.keys()):
                        possive_alt_key.append(temp_possive_alt_key[i])
            
            elif splice_type == "Acceptor-": #o<--
                ref_motif = ref.fetch(sj_key_chr,start-1-1,start+5).upper()

                temp_possive_alt_key = {}
                for i in range(0,len(ref_motif)):
                    if ref_motif[i] in dict3m[i]:
                        possive_alt_motif.append(".")
                    else:
                        if dict3m[i] == "AG":
                            possive_alt_motif.append("R")
                        elif dict3m[i] == "CT":
                            possive_alt_motif.append("Y")
                        else:
                            possive_alt_motif.append(dict3m[i])

                        temp_possive_alt_key[i] = "%d:%s:%s" % (start+i-1, ref_motif[i], dict3m[i])

                if possive_alt_motif[1] == "C" or possive_alt_motif[2] == "T":
                    is_gtag_creation = True
                    if possive_alt_motif[1] == "C":
                        possive_alt_key.append(temp_possive_alt_key[1])

                    if possive_alt_motif[2] == "T":
                        possive_alt_key.append(temp_possive_alt_key[2])
                    
                elif list(set(possive_alt_motif)) == ["."]:
                    possive_alt_key = ['-']
                else:
                    for i in sorted(temp_possive_alt_key.keys()):
                        possive_alt_key.append(temp_possive_alt_key[i])

            else:
                continue

            csvobj["Ref_motif"] = ref_motif
            csvobj["Possive_alt_motif"] = ''.join(possive_alt_motif)
            csvobj["Possive_alt_key"] = ','.join(possive_alt_key)
            csvobj["Is_GTAG_creation"] = is_gtag_creation
            csvwriter.writerow(csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reference = sys.argv[3]
    juncmut_mutpre(input_file, output_file, reference)
