#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 2 2022

@author: Naoko Iida

"""
import csv

def annot_cgd(input_file, output_file, cgd_file):
    
    #make dictionary
    cgd_dict = {}
    with open(cgd_file, 'r') as gin:
        csvreader = csv.DictReader(gin, delimiter='\t')
        for csvobj in csvreader:
            if csvobj["REFERENCES"] != "N/A":
                cgd_dict[csvobj["#GENE"]] = {"CONDITION": csvobj["CONDITION"], "INHERITANCE": csvobj["INHERITANCE"], "AGE GROUP": csvobj["AGE GROUP"]}

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', quotechar='"', fieldnames=csvreader.fieldnames + ["CGD_condition", "CGD_inheritance", "CGD_age"])
        csvwriter.writeheader()

        for csvobj in csvreader:
            condition_list = []
            inheritance_list = []
            age_list = []
            for gene in csvobj["Gene"].split(','):
                if gene in cgd_dict:
                    condition_list.append(cgd_dict[gene]["CONDITION"])
                    inheritance_list.append(cgd_dict[gene]["INHERITANCE"])
                    age_list.append(cgd_dict[gene]["AGE GROUP"])

            if len(condition_list) == 0:
                condition_list.append("---")
                inheritance_list.append("---")
                age_list.append("---")

            csvobj["CGD_condition"] = ','.join(list(set(condition_list)))
            csvobj["CGD_inheritance"] = ','.join(list(set(inheritance_list)))
            csvobj["CGD_age"] = ','.join(list(set(age_list)))

            csvwriter.writerow(csvobj)

if __name__== "__main__":
    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    cgd_file = sys.argv[3]
    
    annot_cgd(input_file, output_file, cgd_file)
