#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_mutpre(input_file, output_file, reference):
    import csv
    import pysam
    
    ref = pysam.FastaFile(reference) # 0-coordinate
    
    trans = str.maketrans('ATGCatgc', 'TACGTAGC')
    
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

            if splice_type == "Donor+": #o-->
                #region -2|0~5 
                bases = ref.fetch(sj_key_chr, start-1-2, start+5).upper()

                mlist = []
                mposi = []
                noncamposi = []
                for i in range(0,len(bases)): 
                    if bases[i] in dict5p[i]:
                        mlist.append(".")
                        mposi.append("%d:%s:%s" % (start+i-2, bases[i], dict5p[i]))
                    else:
                        mlist.append(dict5p[i])
                        mposi.append("%d:%s:%s" % (start+i-2, bases[i], dict5p[i]))
                        noncamposi.append("%d:%s:%s" % (start+i-2, bases[i], dict5p[i]))

                mlist[4] = (mlist[4].replace('AG', 'R'))
                if mlist[2] == "G" or mlist[3] == "T":
                    if mlist[2] == "G" and mlist[3] == "T":
                        mposi_rec = "%s,%s" % (mposi[2], mposi[3])
                        allele = "|%s%s>GT*d" % (bases[2], bases[3])
                    elif mlist[2] == "G" and mlist[3] == ".":
                        mposi_rec = mposi[2]
                        allele = "|%s%s>GT" % (bases[2], bases[3])
                    elif mlist[2] == "." and mlist[3] == "T":
                        mposi_rec = mposi[3]
                        allele = "|%s%s>GT" % (bases[2], bases[3])
                elif len(set(mlist)) == 1:
                    allele = '........'
                    mposi_rec = '-'
                else:
                    allele = 'non-GT/AG'
                    mposi_rec = ','.join(noncamposi)
            
            elif splice_type == "Donor-": #<--o
                bases = ref.fetch(sj_key_chr, end-1-5, end+2).upper()

                mlist = []
                mposi = []
                noncamposi = []
                for i in range(0,len(bases)):
                    if bases[i] in dict5m[i]:
                        mlist.append(".")
                        mposi.append("%d:%s:%s" % (end+i-5, bases[i], dict5m[i]))
                    else:
                        mlist.append(dict5m[i])
                        mposi.append("%d:%s:%s" % (end+i-5, bases[i], dict5m[i]))
                        noncamposi.append("%d:%s:%s" % (end+i-5, bases[i], dict5m[i]))

                mlist[3] = (mlist[3].replace('CT', 'Y'))
                if mlist[4] == "A" or mlist[5] == "C":
                    if mlist[4] == "A" and mlist[5] == "C":
                        mposi_rec = "%s,%s" % (mposi[4], mposi[5])
                        allele = "|%s%s>GT*d" % (bases[5].translate(trans), bases[4].translate(trans))
                    elif mlist[4] == "A" and mlist[5] == ".":
                        mposi_rec = mposi[4]
                        allele = "|%s%s>GT" % (bases[5].translate(trans), bases[4].translate(trans))
                    elif mlist[4] == "." and mlist[5] == "C":
                        mposi_rec = mposi[5] 
                        allele = "|%s%s>GT" % (bases[5].translate(trans), bases[4].translate(trans))
                elif len(set(mlist)) == 1:
                    allele = '........'
                    mposi_rec = '-'
                else:
                    allele = 'non-GT/AG' 
                    mposi_rec = ','.join(noncamposi)
            
            elif splice_type == "Acceptor+": #-->o
                #region -5|1
                bases = ref.fetch(sj_key_chr, end-1-5, end+1).upper()
                mlist = []
                mposi = []
                noncamposi = []
                for i in range(0,len(bases)):
                    if bases[i] in dict3p[i]:
                        mlist.append(".")
                        mposi.append("%d:%s:%s" % (end+i-5, bases[i], dict3p[i]))
                    else:
                        mlist.append(dict3p[i])
                        mposi.append("%d:%s:%s" % (end+i-5, bases[i], dict3p[i]))
                        noncamposi.append("%d:%s:%s" % (end+i-5, bases[i], dict3p[i]))

                mlist[0] = (mlist[0].replace('CT', 'Y'))
                mlist[1] = (mlist[1].replace('CT', 'Y'))
                mlist[3] = (mlist[3].replace('CT', 'Y'))
                mlist[6] = (mlist[6].replace('AG', 'R'))

                if mlist[4] == "A" or mlist[5] == "G":
                    if mlist[4] == "A" and mlist[5] == "G":
                        mposi_rec = "%s,%s" % (mposi[4], mposi[5])
                        allele =  "%s%s|>AG*d" % (bases[4], bases[5])
                    elif mlist[4] == "A" and mlist[5] == ".":
                        mposi_rec = mposi[4]
                        allele =  "%s%s|>AG" % (bases[4], bases[5])
                    elif mlist[4] == "." and mlist[5] == "G":
                        mposi_rec = mposi[5]
                        allele =  "%s%s|>AG" % (bases[4], bases[5])
                    
                elif len(set(mlist)) == 1:
                    allele = '.......'
                    mposi_rec = '-'
                else:
                    allele = 'non-GT/AG' 
                    mposi_rec = ','.join(noncamposi)
            
            elif splice_type == "Acceptor-": #o<--
                bases = ref.fetch(sj_key_chr,start-1-1,start+5).upper()
                mlist = []
                mposi = []
                noncamposi = []
                for i in range(0,len(bases)):
                    if bases[i] in dict3m[i]:
                        mlist.append(".")
                        mposi.append("%d:%s:%s" % (start+i-1, bases[i], dict3m[i]))
                    else:
                        mlist.append(dict3m[i])
                        mposi.append("%d:%s:%s" % (start+i-1, bases[i], dict3m[i]))
                        noncamposi.append("%d:%s:%s" % (start+i-1, bases[i], dict3m[i]))

                mlist[0] = (mlist[0].replace('CT', 'Y'))
                mlist[3] = (mlist[3].replace('AG', 'R'))
                mlist[5] = (mlist[5].replace('AG', 'R'))
                mlist[6] = (mlist[6].replace('AG', 'R'))

                if mlist[1] == "C" or mlist[2] == "T":
                    if mlist[1] == "C" and mlist[2] == "T":
                        mposi_rec = "%s,%s" % (mposi[1], mposi[2])
                        allele =  "%s%s|>AG*d" % (bases[2].translate(trans), bases[1].translate(trans))
                    elif mlist[1] == "C" and mlist[2] == ".":
                        mposi_rec = mposi[1]
                        allele =  "%s%s|>AG" % (bases[2].translate(trans), bases[1].translate(trans))
                    elif mlist[1] == "." and mlist[2] == "T":
                        mposi_rec = mposi[2]
                        allele =  "%s%s|>AG" % (bases[2].translate(trans), bases[1].translate(trans))
                    
                elif len(set(mlist)) == 1:
                    allele = '.......'
                    mposi_rec = '-'
                else:
                    allele = 'non-GT/AG' 
                    mposi_rec = ','.join(noncamposi)

            else:
                continue

            csvobj["Ref_motif"] = bases
            csvobj["Possive_alt_motif"] = ''.join(mlist)
            csvobj["Possive_alt_key"] = mposi_rec
            csvobj["Is_GTAG_creation"] = allele
            csvwriter.writerow(csvobj)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reference = sys.argv[3]
    juncmut_mutpre(input_file, output_file, reference)
