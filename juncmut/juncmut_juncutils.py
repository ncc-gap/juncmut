#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#! /usr/bin/env python

import pysam 
import sys
import csv
import subprocess
import shutil
import os
import gzip

# for python2 and 3 compatibility
try:
    import itertools.izip as zip
except ImportError:
    pass

import gzip, subprocess 

def gtf_to_bed(output_gene_file, output_exon_file, gencode_gene_file):
    if gencode_gene_file.endswith("gtf.gz"):
        key_value_split = " "
    elif gencode_gene_file.endswith("gff3.gz"):
        key_value_split = "="
    else:
        raise Exception("unexcepted gencode gene file")

    with gzip.open(gencode_gene_file, 'rt') as hin, \
        open(output_gene_file + ".unsorted.tmp", 'w') as hout_gene, \
        open(output_exon_file + ".unsorted.tmp", 'w') as hout_exon:

        csvwriter_gene = csv.DictWriter(hout_gene, delimiter='\t', lineterminator='\n', fieldnames=[
            "Chr", "Gene_start", "Gene_end", "Gene_name", "Score", "Strand", "Symbol", "MANE"
        ])

        csvwriter_exon = csv.DictWriter(hout_exon, delimiter='\t', lineterminator='\n', fieldnames=[
            "Chr", "Exon_start", "Exon_end", "Gene_name", "Exon_num", "Strand"
        ])

        for line in hin:
            if line.startswith("#"):
                continue
            F = line.rstrip('\n').split('\t')
            feature_type = F[2]
            if feature_type != "transcript" and feature_type != "exon":
                continue

            chr = F[0]
            start = int(F[3]) - 1
            end = F[4]
            strand = F[6]

            add_info = {"MANE": "---", }
            for item in F[8].rstrip(";").split(";"):
                (key, value) = item.strip(" ").replace('"', '').split(key_value_split)
                if key in ["transcript_id", "gene_name", "exon_number"]:
                    add_info[key] = value
                elif key == "tag" and "MANE" in value:
                    add_info["MANE"] = value

            if feature_type == "transcript":
                csvwriter_gene.writerow({
                    "Chr": chr,
                    "Gene_start": start,
                    "Gene_end": end,
                    "Gene_name": add_info["transcript_id"],
                    "Score": 0, 
                    "Strand": strand,
                    "Symbol": add_info["gene_name"],
                    "MANE": add_info["MANE"],
                })

            else:
                csvwriter_exon.writerow({
                    "Chr": chr,
                    "Exon_start": start,
                    "Exon_end": end,
                    "Gene_name": add_info["transcript_id"],
                    "Exon_num": add_info["exon_number"], 
                    "Strand": strand
                })

    with open(output_gene_file + ".sorted.tmp", 'w') as hout_gene:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_gene_file + ".unsorted.tmp"], stdout = hout_gene)

    with open(output_gene_file, 'w') as hout_gene:
        subprocess.check_call(["bgzip", "-f", "-c", output_gene_file + ".sorted.tmp"], stdout = hout_gene)
    
    subprocess.check_call(["tabix", "-p", "bed", output_gene_file])
    subprocess.check_call(["rm", "-rf", output_gene_file + ".unsorted.tmp"])
    subprocess.check_call(["rm", "-rf", output_gene_file + ".sorted.tmp"])

    with open(output_exon_file + ".sorted.tmp", 'w') as hout_exon:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_exon_file + ".unsorted.tmp"], stdout = hout_exon)

    with open(output_exon_file, 'w') as hout_exon:
        subprocess.check_call(["bgzip", "-f", "-c", output_exon_file + ".sorted.tmp"], stdout = hout_exon)

    subprocess.check_call(["tabix", "-p", "bed", output_exon_file])
    subprocess.check_call(["rm", "-rf", output_exon_file + ".unsorted.tmp"])
    subprocess.check_call(["rm", "-rf", output_exon_file + ".sorted.tmp"])


def proc_star_junction(input_file, output_file, control_file, read_num_thres, overhang_thres, remove_annotated, convert_map_splice2):

    is_control = True if control_file is not None else False

    control_db = {}
    if is_control:
        with gzip.open(control_file, 'rt') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = F[0] + '\t' + F[1] + '\t' + F[2]
                control_db[key] = 1


    if read_num_thres is None: read_num_thres = 0
    if overhang_thres is None: overhang_thres = 0
    if remove_annotated is None: remove_annotated = False
    
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2]
            if remove_annotated == True and F[5] != "0": continue
            if int(F[6]) < read_num_thres: continue
            if int(F[8]) < overhang_thres: continue


            if key in control_db: continue
            """
            ##########
            # remove control files
            if is_control:
                tabixErrorFlag = 0
                try:
                    records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    # tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                control_flag = 0;
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                            control_flag = 1

                if control_flag == 1: continue
            ##########
            """

            if convert_map_splice2:
                # convert to map-splice2 coordinate
                F[1] = str(int(F[1]) - 1)
                F[2] = str(int(F[2]) + 1)

            print('\t'.join(F), file = hout)

    hout.close()

def juncutils_annotate(input_file, output_file, gencode_gene_file, junction_margin, exon_margin, debug=False):
    """
        The purpose of this script is to classify splicing changes
        mainly by comparing the two breakpoints with the exon-intorn junction of genes
        within the database.
        Also, we generate the sequence arrond the breakpoints, which will be helpful
        for checking the authenticity of the splicing and evaluating the relationships
        with the somatic mutations.

        here is the classification categories:
        1. known (The splicing pattern is included in the database)
        (the start and end breakpoints are next exon-intron junctions of the same gene) 
        2. exon skipping 
        (the start and end breakpoints are exon-intron junctions of the same gene,
         but not the next ones)
        3. splice-site slip
        (one of the two breakpoints is an exon-intron junction and the other is within the 30bp exon of the same gene)
        4. pseudo-exon inclusion
        (one of the two break points is an exon-intron junction and the other is located in the same gene, but more than 30bp from exons of the gene)
        5. other
        (neighter of the two breakpoins are exon-intron junction, but located in the same gene)
        6. chimeric (spliced)
        7. chimeric (un-spliced)


        The algorithm for the annotation is as follows
        1. for both breakpoints, list up the exon-intron junctions matching to the breakpoints
        2. for both breakpoints, list up the exons within 30bp from the breakpoints
        3. for both breakpoints, list up the genes matching to the breakpoints
        4. summarize the above results and induce the annotation from them
        5. get the sequence arround the breakpoints.

    """
    
    gtf_to_bed(output_file + ".tmp.refGene.bed.gz", output_file + ".tmp.refExon.bed.gz", gencode_gene_file)

    with open(output_file + ".tmp1.junc1.gene.bed", 'w') as hout_g1, open(output_file + ".tmp1.junc2.gene.bed", 'w') as hout_g2, \
      open(output_file + ".tmp1.junc1.exon.bed", 'w') as hout_e1, open(output_file + ".tmp1.junc2.exon.bed", 'w') as hout_e2:
        with open(input_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                chr_name, sj_start, sj_end = F[0], int(F[1]) - 1, int(F[2]) + 1
                sj_id = ','.join([F[0], F[1], F[2]])
                
                print(f'{chr_name}\t{max(0, sj_start - 1)}\t{sj_start}\t{sj_id}', file = hout_g1)
                print(f'{chr_name}\t{max(0, sj_end - 1)}\t{sj_end}\t{sj_id}', file = hout_g2)
                print(f'{chr_name}\t{max(0, sj_start - exon_margin - 1)}\t{sj_start + exon_margin}\t{sj_id}', file = hout_e1)
                print(f'{chr_name}\t{max(0, sj_end - exon_margin - 1)}\t{sj_end + exon_margin}\t{sj_id}', file = hout_e2) 

    with open(output_file + ".tmp2.junc1.gene.bed", 'w') as hout_g1:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc1.gene.bed", "-b", output_file + ".tmp.refGene.bed.gz", "-loj"], stdout = hout_g1)

    with open(output_file + ".tmp2.junc2.gene.bed", 'w') as hout_g2:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc2.gene.bed", "-b", output_file + ".tmp.refGene.bed.gz", "-loj"], stdout = hout_g2)

    with open(output_file + ".tmp2.junc1.exon.bed", 'w') as hout_e1:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc1.exon.bed", "-b", output_file + ".tmp.refExon.bed.gz", "-loj"], stdout = hout_e1)

    with open(output_file + ".tmp2.junc2.exon.bed", 'w') as hout_e2:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc2.exon.bed", "-b", output_file + ".tmp.refExon.bed.gz", "-loj"], stdout = hout_e2)

    if not debug:
        subprocess.check_call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc1.gene.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc2.gene.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc1.exon.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc2.exon.bed"])

       
    with open(output_file + ".tmp2.junc1.gene.bed", 'r') as hin, open(output_file + ".tmp3.junc1.gene.bed", 'w') as hout:
        tmp_id, tmp_gene, tmp_symbol, tmp_mane = "", [], [], []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] != tmp_id:
                if tmp_id != "": print(tmp_id + '\t' + ','.join(tmp_gene) + '\t' + ','.join(tmp_symbol) + '\t' + ','.join(tmp_mane), file = hout)
                tmp_id, tmp_gene, tmp_symbol, tmp_mane = F[3], [], [], []
            if F[7] != ".":
                if not F[7] in tmp_gene:
                    tmp_gene.append(F[7])
                    tmp_symbol.append(F[10])
                    tmp_mane.append(F[11])
            else:
                tmp_gene.append("---")
                tmp_symbol.append("---")
                tmp_mane.append("---")

        if tmp_id != "": print(tmp_id + '\t' + ','.join(tmp_gene) + '\t' + ','.join(tmp_symbol) + '\t' + ','.join(tmp_mane), file = hout)

    with open(output_file + ".tmp2.junc2.gene.bed", 'r') as hin, open(output_file + ".tmp3.junc2.gene.bed", 'w') as hout:
        tmp_id, tmp_gene, tmp_symbol, tmp_mane = "", [], [], []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] != tmp_id:
                if tmp_id != "": print(tmp_id + '\t' + ','.join(tmp_gene) + '\t' + ','.join(tmp_symbol) + '\t' + ','.join(tmp_mane), file = hout)
                tmp_id, tmp_gene, tmp_symbol, tmp_mane = F[3], [], [], []
            if F[7] != ".":
                if not F[7] in tmp_gene:
                    tmp_gene.append(F[7])
                    tmp_symbol.append(F[10])
                    tmp_mane.append(F[11])
            else:
                tmp_gene.append("---")
                tmp_symbol.append("---")
                tmp_mane.append("---")

        if tmp_id != "": print(tmp_id + '\t' + ','.join(tmp_gene) + '\t' + ','.join(tmp_symbol) + '\t' + ','.join(tmp_mane), file = hout)

    with open(output_file + ".tmp2.junc1.exon.bed", 'r') as hin, open(output_file + ".tmp3.junc1.exon.bed", 'w') as hout:
        tmp_id, tmp_gene, tmp_exon_num, tmp_edge, tmp_offset = "", [], [], [], []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            FF = F[3].split(',')
            chr_name, sj_start, sj_end = FF[0], int(FF[1]) - 1, int(FF[2]) + 1

            if F[3] != tmp_id:
                if tmp_id != "":
                    if not len(tmp_gene) == len(tmp_exon_num) == len(tmp_edge) == len(tmp_offset):
                        print("Inconsistency for the format in creating exon information files", file = sys.stderr)
                        sys.exit(1)
                    print('\t'.join([tmp_id, ','.join(tmp_gene), ','.join(tmp_exon_num), ','.join(tmp_edge), ','.join(tmp_offset)]), file = hout)
                tmp_id, tmp_gene, tmp_exon_num, tmp_edge, tmp_offset = F[3], [], [], [], []

            if F[7] != ".":
                tmp_gene.append(F[7])
                tmp_exon_num.append(F[8])
                if abs(sj_start - int(F[6])) < junction_margin:
                    tmp_offset.append(str(sj_start - int(F[6])))
                    if F[9] == '+': tmp_edge.append('e')
                    if F[9] == '-': tmp_edge.append('s')
                else:
                    tmp_edge.append("---"), tmp_offset.append("---")
            else:
                tmp_gene.append("---"), tmp_exon_num.append("---"), tmp_edge.append("---"), tmp_offset.append("---")

        if tmp_id != "":
            if not len(tmp_gene) == len(tmp_exon_num) == len(tmp_edge) == len(tmp_offset):
                print("Inconsistency for the format in creating exon information files", file = sys.stderr)
                sys.exit(1)
            print('\t'.join([tmp_id, ','.join(tmp_gene), ','.join(tmp_exon_num), ','.join(tmp_edge), ','.join(tmp_offset)]), file = hout)

    with open(output_file + ".tmp2.junc2.exon.bed", 'r') as hin, open(output_file + ".tmp3.junc2.exon.bed", 'w') as hout:
        tmp_id, tmp_gene, tmp_exon_num, tmp_edge, tmp_offset = "", [], [], [], []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            FF = F[3].split(',')
            chr_name, sj_start, sj_end = FF[0], int(FF[1]) - 1, int(FF[2]) + 1

            if F[3] != tmp_id:
                if tmp_id != "":
                    if not len(tmp_gene) == len(tmp_exon_num) == len(tmp_edge) == len(tmp_offset):
                        print("Inconsistency for the format in creating exon information files", file = sys.stderr)
                        sys.exit(1)
                    print('\t'.join([tmp_id, ','.join(tmp_gene), ','.join(tmp_exon_num), ','.join(tmp_edge), ','.join(tmp_offset)]), file = hout)
                tmp_id, tmp_gene, tmp_exon_num, tmp_edge, tmp_offset = F[3], [], [], [], []

            if F[7] != ".":
                tmp_gene.append(F[7])
                tmp_exon_num.append(F[8])
                if abs(sj_end - 1 - int(F[5])) < junction_margin:
                    tmp_offset.append(str(sj_end - 1 - int(F[5])))
                    if F[9] == '+': tmp_edge.append('s')
                    if F[9] == '-': tmp_edge.append('e')
                else:
                    tmp_edge.append("---"), tmp_offset.append("---")
            else:
                tmp_gene.append("---"), tmp_exon_num.append("---"), tmp_edge.append("---"), tmp_offset.append("---")

        if tmp_id != "":
            if not len(tmp_gene) == len(tmp_exon_num) == len(tmp_edge) == len(tmp_offset):
                print("Inconsistency for the format in creating exon information files", file = sys.stderr)
                sys.exit(1)
            print('\t'.join([tmp_id, ','.join(tmp_gene), ','.join(tmp_exon_num), ','.join(tmp_edge), ','.join(tmp_offset)]), file = hout)

    if not debug:
        subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc1.gene.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc2.gene.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc1.exon.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc2.exon.bed"])

    ##########


    with open(input_file, 'r') as hin:
        input_header = hin.readline().rstrip('\n').split('\t')

    with open(input_file, 'r') as hin, open(output_file + ".tmp3.junc1.gene.bed", 'r') as hin_g1, open(output_file + ".tmp3.junc2.gene.bed", 'r') as hin_g2, \
      open(output_file + ".tmp3.junc1.exon.bed", 'r') as hin_e1, open(output_file + ".tmp3.junc2.exon.bed", 'r') as hin_e2, \
      open(output_file, 'w') as hout:

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=["SJ_" + str(i) for i in range(1, len(input_header) + 1)] + [
#            "Splicing_Class", "Is_Inframe", "Gene_1", "Symbol_1", "Exon_Num_1", "Is_Boundary_1", "Offset_1", "Gene_2", "Symbol_2", "Exon_Num_2", "Is_Boundary_2", "Offset_2"
            "Splicing_Class", "Gene_1", "Symbol_1", "MANE_1", "Is_Boundary_1", "Offset_1", "Gene_2", "Symbol_2", "MANE_2", "Is_Boundary_2", "Offset_2"
        ])
        csvwriter.writeheader()

        for line, line_g1, line_g2, line_e1, line_e2 in zip(hin, hin_g1, hin_g2, hin_e1, hin_e2):

            F = line.rstrip('\n').split('\t')
            F_g1, F_g2, F_e1, F_e2 = line_g1.rstrip('\n').split('\t'), line_g2.rstrip('\n').split('\t'), \
                                     line_e1.rstrip('\n').split('\t'), line_e2.rstrip('\n').split('\t')

            chr_name, sj_start, sj_end = F[0], int(F[1]) - 1, int(F[1]) + 1
            sj_id = ','.join([F[0], F[1], F[2]])
            ##########

            # key check
            if not sj_id == F_g1[0] == F_g2[0] == F_e1[0] == F_e2[0]:
                print("Inconsistency of splicing junction keys in the information files", file = sys.stderr)
                print(sj_id, F_g1[0], F_g2[0], F_e1[0], F_e2[0])
                sys.exit(1)

            gene1 = F_g1[1].split(',') if F_g1[1] != "---" else []
            gene2 = F_g2[1].split(',') if F_g2[1] != "---" else []
            symbol1 = F_g1[2].split(',') if F_g1[2] != "---" else []
            symbol2 = F_g2[2].split(',') if F_g2[2] != "---" else []
            mane1 = F_g1[3].split(',') if F_g1[2] != "---" else []
            mane2 = F_g2[3].split(',') if F_g2[2] != "---" else []

            exon1, junction1, offset1 = {}, {}, {}
            for tmp_gene, tmp_exon_num, tmp_edge, tmp_offset in zip(F_e1[1].split(','), F_e1[2].split(','), F_e1[3].split(','), F_e1[4].split(',')):
                if tmp_gene == "---": continue
                exon1[tmp_gene] = int(tmp_exon_num)
                if tmp_edge != "---": junction1[tmp_gene] = tmp_edge
                if tmp_offset != "---": offset1[tmp_gene] = tmp_offset
                 
            exon2, junction2, offset2 = {}, {}, {}
            for tmp_gene, tmp_exon_num, tmp_edge, tmp_offset in zip(F_e2[1].split(','), F_e2[2].split(','), F_e2[3].split(','), F_e2[4].split(',')):
                if tmp_gene == "---": continue
                exon2[tmp_gene] = int(tmp_exon_num)
                if tmp_edge != "---": junction2[tmp_gene] = tmp_edge
                if tmp_offset != "---": offset2[tmp_gene] = tmp_offset

            spliceClass = ""
            checkGenes = gene1 + gene2

            ##########
            # check for know junction or exon skip
            spliceClass = ""
            for gene_index, gene in enumerate(checkGenes):
                if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
                    if (junction1[gene] == "e" and junction2[gene] == "s" and exon2[gene] - exon1[gene] > 0) or \
                       (junction2[gene] == "e" and junction1[gene] == "s" and exon1[gene] - exon2[gene] > 0):
                        spliceClass = "Known or Exon skipping"
                        break

            if spliceClass != "":
                #out_csvobj = {
                #    "Splicing_Class": "Known or Exon skipping", 
                #    "Gene_1": '---', 
                #    "Symbol_1": '---', 
                #    "Is_Boundary_1": '---', 
                #    "Offset_1": '---', 
                #    "Gene_2": '---', 
                #    "Symbol_2": '---', 
                #    "Is_Boundary_2": '---', 
                #    "Offset_2": '---'
                #}
                #for i in range(0, len(F)):
                #    out_csvobj["SJ_" + str(i+1)] = F[i]
                #csvwriter.writerow(out_csvobj)
                continue

            ##########
            # check for alternative-3'-splice-site or alternative-5'-splice-site
            passGene_acceptor = []
            passGene_donor = []
            for gene_index, gene in enumerate(checkGenes):
                if gene in gene1 and gene in gene2:
                    if (gene in junction1 and junction1[gene] == "e" and gene not in junction2) or \
                       (gene in junction2 and junction2[gene] == "e" and gene not in junction1):
                        passGene_acceptor.append(gene_index)

                    if (gene in junction1 and junction1[gene] == "s" and gene not in junction2) or \
                       (gene in junction2 and junction2[gene] == "s" and gene not in junction1):
                        passGene_donor.append(gene_index)

            passGene = []
            if len(passGene_acceptor) > 0:
                spliceClass = "Alternative 3'SS"
                passGene.extend(passGene_acceptor)
            elif len(passGene_donor) > 0:
                spliceClass = "Alternative 5'SS"
                passGene.extend(passGene_donor)
            else:
                #out_csvobj = {
                #    "Splicing_Class": "Other", 
                #    "Gene_1": '---', 
                #    "Symbol_1": '---', 
                #    "Is_Boundary_1": '---', 
                #    "Offset_1": '---', 
                #    "Gene_2": '---', 
                #    "Symbol_2": '---', 
                #    "Is_Boundary_2": '---', 
                #    "Offset_2": '---'
                #}
                #for i in range(0, len(F)):
                #    out_csvobj["SJ_" + str(i+1)] = F[i]
                #csvwriter.writerow(out_csvobj)
                continue
            

            # summarize the exon and junction information for display
            symbols = symbol1 + symbol2
            manes = mane1 + mane2

            geneInfo1 = []
            symbolInfo1 = []
            maneInfo1 = []
            junctionInfo1 = []
            offsetInfo1 = []
            if len(gene1) > 0:
                for g1 in sorted(gene1):
                    g1_index = checkGenes.index(g1)
                    if g1_index not in passGene: continue 
                    geneInfo1.append(g1)
                    symbolInfo1.append(symbols[g1_index])
                    maneInfo1.append(manes[g1_index])

                    if g1 in junction1:
                        junctionInfo1.append(junction1[g1])
                    else:
                        junctionInfo1.append("*")

                    if g1 in offset1:
                        offsetInfo1.append(str(offset1[g1]))
                    else:
                        offsetInfo1.append("*")

            if len(geneInfo1) == 0: 
                geneInfo1.append("---")
                symbolInfo1.append("---")
                maneInfo1.append("---")
                junctionInfo1.append("---")
                offsetInfo1.append("---")

            geneInfo2 = []
            symbolInfo2 = []
            maneInfo2 = []
            junctionInfo2 = []
            offsetInfo2 = []
            if len(gene2) > 0:
                for g2 in sorted(gene2):
                    g2_index = checkGenes.index(g2)
                    if g2_index not in passGene: continue
                    geneInfo2.append(g2)
                    symbolInfo2.append(symbols[g2_index])
                    maneInfo2.append(manes[g2_index])

                    if g2 in junction2:
                        junctionInfo2.append(junction2[g2])
                    else:
                        junctionInfo2.append("*")

                    if g2 in offset2:
                        offsetInfo2.append(str(offset2[g2]))
                    else:
                        offsetInfo2.append("*")

            if len(geneInfo2) == 0:
                geneInfo2.append("---")
                symbolInfo2.append("---")
                maneInfo2.append("---")
                junctionInfo2.append("---")
                offsetInfo2.append("---")

         
            out_csvobj = {
                "Splicing_Class": spliceClass, 
                #"Is_Inframe": in_frame, 
                "Gene_1": ';'.join(geneInfo1), 
                "Symbol_1": ';'.join(symbolInfo1), 
                "MANE_1": ';'.join(maneInfo1), 
                #"Exon_Num_1": ';'.join(exonInfo1), 
                "Is_Boundary_1": ';'.join(junctionInfo1), 
                "Offset_1": ';'.join(offsetInfo1), 
                "Gene_2": ';'.join(geneInfo2), 
                "Symbol_2": ';'.join(symbolInfo2), 
                "MANE_2": ';'.join(maneInfo2), 
                #"Exon_Num_2": ';'.join(exonInfo2), 
                "Is_Boundary_2": ';'.join(junctionInfo2), 
                "Offset_2": ';'.join(offsetInfo2)
            }
            for i in range(0, len(F)):
                out_csvobj["SJ_" + str(i+1)] = F[i]
            csvwriter.writerow(out_csvobj)

    if not debug:
        subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc1.gene.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc2.gene.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc1.exon.bed"])
        subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc2.exon.bed"])

def juncmut_juncutils(input_file, output_file, control_list, genecode_gene_file, read_num_thres, junction_margin, exon_margin, debug):

    target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

    tmpfile1 = output_file + ".tmp1"
    with open(input_file, 'r') as hin, open(tmpfile1, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] in target_rnames + ["chr" + x for x in target_rnames]:
                print('\t'.join(F), file = hout)

    tmpfile_list = []
    tmpfile_list.append(tmpfile1)

    tmpfile2 = output_file + ".tmp2"
    if not control_list:
        proc_star_junction(tmpfile1, tmpfile2, None, read_num_thres, 10, False, False)
        tmpfile_list.append(tmpfile2)
    else:
        cur_infile = tmpfile1
        for n,control in enumerate(control_list):
            cur_outfile = "%s_%d" % (tmpfile2, n)
            proc_star_junction(cur_infile, cur_outfile, control, read_num_thres, 10, False, False)
            tmpfile_list.append(cur_outfile)
            cur_infile = cur_outfile

        shutil.copy(cur_infile, tmpfile2)
        tmpfile_list.append(tmpfile2)

    juncutils_annotate(tmpfile2, output_file, genecode_gene_file, junction_margin, exon_margin, debug)

    if not debug:
        for tfile in tmpfile_list:
            os.remove(tfile)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    control_list = sys.argv[3].split(",")
    genecode_gene_file = sys.argv[4]

    juncmut_juncutils(input_file, output_file, control_list, genecode_gene_file, 1, 3, 30, True)
