#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#! /usr/bin/env python

import pysam 
import sys, subprocess

# for python2 and 3 compatibility
try:
    import itertools.izip as zip
except ImportError:
    pass

import annot_utils.gene, annot_utils.exon, annot_utils.coding

def juncutils_annotate(input_file, output_file, ucsc_gene_file, junction_margin, exon_margin, gene_model = "refseq", genome_id = "hg19", is_grc = True):
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

    annot_utils.gene.make_gene_info(output_file + ".tmp.refGene.bed.gz", ucsc_gene_file, gene_model, genome_id, is_grc, True)
    annot_utils.exon.make_exon_info(output_file + ".tmp.refExon.bed.gz", ucsc_gene_file, gene_model, genome_id, is_grc, True)
    annot_utils.coding.make_coding_info(output_file + ".tmp.refCoding.bed.gz", ucsc_gene_file, gene_model, genome_id, is_grc, True)
  
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
                """
                print(chr_name + '\t' + str(sj_start - 1) + '\t' + str(sj_start) + '\t' + sj_id, file = hout_g1)
                print(chr_name + '\t' + str(sj_end - 1) + '\t' + str(sj_end) + '\t' + sj_id, file = hout_g2)
                print(chr_name + '\t' + str(sj_start - exon_margin - 1) + '\t' + str(sj_start + exon_margin) + '\t' + sj_id, file = hout_e1)
                print(chr_name + '\t' + str(sj_end - exon_margin - 1) + '\t' + str(sj_end + exon_margin) + '\t' + sj_id, file = hout_e2)
                """

    with open(output_file + ".tmp2.junc1.gene.bed", 'w') as hout_g1:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc1.gene.bed", "-b", output_file + ".tmp.refGene.bed.gz", "-loj"], stdout = hout_g1)

    with open(output_file + ".tmp2.junc2.gene.bed", 'w') as hout_g2:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc2.gene.bed", "-b", output_file + ".tmp.refGene.bed.gz", "-loj"], stdout = hout_g2)

    with open(output_file + ".tmp2.junc1.exon.bed", 'w') as hout_e1:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc1.exon.bed", "-b", output_file + ".tmp.refExon.bed.gz", "-loj"], stdout = hout_e1)

    with open(output_file + ".tmp2.junc2.exon.bed", 'w') as hout_e2:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp1.junc2.exon.bed", "-b", output_file + ".tmp.refExon.bed.gz", "-loj"], stdout = hout_e2)

    subprocess.check_call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc1.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc2.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc1.exon.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp1.junc2.exon.bed"])

       
    with open(output_file + ".tmp2.junc1.gene.bed", 'r') as hin, open(output_file + ".tmp3.junc1.gene.bed", 'w') as hout:
        tmp_id, tmp_gene = "", []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] != tmp_id:
                if tmp_id != "": print(tmp_id + '\t' + ','.join(list(set(tmp_gene))), file = hout)
                tmp_id, tmp_gene = F[3], []
            if F[7] != ".":
                tmp_gene.append(F[7])
            else:
                tmp_gene.append("---")

        if tmp_id != "": print(tmp_id + '\t' + ','.join(list(set(tmp_gene))), file = hout)

    with open(output_file + ".tmp2.junc2.gene.bed", 'r') as hin, open(output_file + ".tmp3.junc2.gene.bed", 'w') as hout:
        tmp_id, tmp_gene = "", []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] != tmp_id:
                if tmp_id != "": print(tmp_id + '\t' + ','.join(list(set(tmp_gene))), file = hout)
                tmp_id, tmp_gene = F[3], []
            if F[7] != ".":
                tmp_gene.append(F[7])
            else:
                tmp_gene.append("---")

        if tmp_id != "": print(tmp_id + '\t' + ','.join(list(set(tmp_gene))), file = hout)

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

    subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc1.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc2.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc1.exon.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp2.junc2.exon.bed"])

    ##########
    # for creating splicing junction to in-frame information table
    with open(output_file + ".tmp3.junc1.gene.bed") as hin_g1, open(output_file + ".tmp3.junc2.gene.bed") as hin_g2, \
      open(output_file + ".tmp3.junc12.gene.bed", 'w') as hout:
        for line_g1, line_g2 in zip(hin_g1, hin_g2):
            F_g1 = line_g1.rstrip('\n').split('\t')
            F_g2 = line_g2.rstrip('\n').split('\t')

            # key check
            if not F_g1[0] == F_g2[0]:
                print("Inconsistency of splicing junction keys in the information files: " + output_file + ".tmp3.junc12.gene.bed", file = sys.stderr)
                print(F_g1[0] + ' ' + F_g2[0])
                sys.exit(1)

            FF = F_g1[0].split(',')
            chr_name, sj_start, sj_end = FF[0], int(FF[1]) - 1, int(FF[2]) + 1

            gene1 = F_g1[1].split(',') if F_g1[1] != "---" else []
            gene2 = F_g2[1].split(',') if F_g2[1] != "---" else []

            common_gene_list = list(set(gene1) & set(gene2))
            common_gene = ','.join(common_gene_list) if len(common_gene_list) > 0 else "---"

            print(f'{chr_name}\t{max(0, sj_start - exon_margin - 1)}\t{sj_end + exon_margin}\t{F_g1[0]}\t{common_gene}', file = hout)

            """
            print(chr_name + '\t' + str(sj_start - exon_margin - 1) + '\t' + str(sj_end + exon_margin) + '\t' + F_g1[0] + '\t' + common_gene, file = hout)
            """

    with open(output_file + ".tmp3.junc12.gene.coding.bed", 'w') as hout:
        subprocess.check_call(["bedtools", "intersect", "-a", output_file + ".tmp3.junc12.gene.bed", "-b", output_file + ".tmp.refCoding.bed.gz", "-loj"], stdout = hout)

    with open(output_file + ".tmp3.junc12.gene.coding.bed") as hin, open(output_file + ".tmp3.junc12.coding_size.bed", 'w') as hout:
        tmp_id, tmp_sj_start, tmp_sj_end, tmp_gene2coding_size = "", 0, 0, {}
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if tmp_id != F[3]:
                if tmp_id != "":
                    if len(tmp_gene2coding_size) == 0: tmp_gene2coding_size["---"] = "---"
                    print(tmp_id + '\t' + ','.join(tmp_gene2coding_size.keys()) + '\t' + \
                      ','.join([str(x) for x in tmp_gene2coding_size.values()]), file = hout)

                tmp_id = F[3]
                _, tmp_sj_start, tmp_sj_end = tmp_id.split(',')
                tmp_sj_start, tmp_sj_end = int(tmp_sj_start) - 1, int(tmp_sj_end) + 1 
                tmp_gene2coding_size = {}

            if F[4] == "---": continue
            if F[9] != "coding": continue
            if F[8] not in F[4].split(','): continue

            exon_start, exon_end, exon_gene = int(F[6]), int(F[7]), F[8]
            if exon_gene not in tmp_gene2coding_size: tmp_gene2coding_size[exon_gene] = 0

            # sj_start overlaps with the current coding region
            if exon_start <= tmp_sj_start + exon_margin and tmp_sj_start - exon_margin <= exon_end:
                # sj_end also overlaps with the current coding region
                if exon_start <= tmp_sj_end + exon_margin and tmp_sj_end - exon_margin <= exon_end:
                    tmp_gene2coding_size[exon_gene] = tmp_gene2coding_size[exon_gene] + tmp_sj_end - tmp_sj_start - 1
                else:
                    tmp_gene2coding_size[exon_gene] = tmp_gene2coding_size[exon_gene] + exon_end - tmp_sj_start
            # sj_end overlaps with the current coding region
            elif exon_start <= tmp_sj_end + exon_margin and tmp_sj_end - exon_margin <= exon_end:
                tmp_gene2coding_size[exon_gene] = tmp_gene2coding_size[exon_gene] + tmp_sj_end - 1 - exon_start 
            else:                
                tmp_gene2coding_size[exon_gene] = tmp_gene2coding_size[exon_gene] + exon_end - exon_start

        if tmp_id != "":
            if len(tmp_gene2coding_size) == 0: tmp_gene2coding_size["---"] = "---"
            print(tmp_id + '\t' + ','.join(tmp_gene2coding_size.keys()) + '\t' + \
                  ','.join([str(x) for x in tmp_gene2coding_size.values()]), file = hout)
 
    subprocess.check_call(["rm", "-rf", output_file + ".tmp.refCoding.bed.gz"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp.refCoding.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc12.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc12.gene.coding.bed"])
    ##########


    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        line = hin.readline()
        F = line.rstrip('\n').split('\t')
        print('\t'.join(["SJ_" + str(i) for i in range(1, len(F) + 1)]) + '\t' + \
              "Splicing_Class" + '\t' + "Is_Inframe" + '\t' + "Gene_1" + '\t' + "Exon_Num_1" + '\t' + \
              "Is_Boundary_1" + '\t' + "Offset_1" + '\t' + "Gene_2" + '\t' + "Exon_Num_2" + '\t' + \
              "Is_Boundary_2" + '\t' + "Offset_2", file = hout)

    with open(input_file, 'r') as hin, open(output_file + ".tmp3.junc1.gene.bed", 'r') as hin_g1, open(output_file + ".tmp3.junc2.gene.bed", 'r') as hin_g2, \
      open(output_file + ".tmp3.junc1.exon.bed", 'r') as hin_e1, open(output_file + ".tmp3.junc2.exon.bed", 'r') as hin_e2, \
      open(output_file + ".tmp3.junc12.coding_size.bed") as hin_c, open(output_file, 'a') as hout:

        for line, line_g1, line_g2, line_e1, line_e2, line_c in zip(hin, hin_g1, hin_g2, hin_e1, hin_e2, hin_c):

            F = line.rstrip('\n').split('\t')
            F_g1, F_g2, F_e1, F_e2, F_c = line_g1.rstrip('\n').split('\t'), line_g2.rstrip('\n').split('\t'), \
                                            line_e1.rstrip('\n').split('\t'), line_e2.rstrip('\n').split('\t'), line_c.rstrip('\n').split('\t')

            chr_name, sj_start, sj_end = F[0], int(F[1]) - 1, int(F[1]) + 1
            sj_id = ','.join([F[0], F[1], F[2]])
            ##########

            # key check
            if not sj_id == F_g1[0] == F_g2[0] == F_e1[0] == F_e2[0] == F_c[0]:
                print("Inconsistency of splicing junction keys in the information files", file = sys.stderr)
                print(sj_id, F_g1[0], F_g2[0], F_e1[0], F_e2[0], F_c[0])
                sys.exit(1)

            gene1 = F_g1[1].split(',') if F_g1[1] != "---" else []
            gene2 = F_g2[1].split(',') if F_g2[1] != "---" else [] 

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


            gene2coding_size = {}
            for tmp_gene, tmp_coding_size in zip(F_c[1].split(','), F_c[2].split(',')):
                if tmp_gene == "---": continue
                gene2coding_size[tmp_gene] = int(tmp_coding_size)
 

            spliceClass = ""
            in_frame = "---"
            checkGenes = list(set(gene1 + gene2))
            ##########
            # check for know junction
            passGene = []
            for gene in checkGenes:
                if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
                    if junction1[gene] == "e" and junction2[gene] == "s" and exon2[gene] - exon1[gene] == 1: passGene.append(gene)
                    if junction2[gene] == "e" and junction1[gene] == "s" and exon1[gene] - exon2[gene] == 1: passGene.append(gene)

            if len(passGene) > 0: spliceClass = "Known"

            ##########
            # check for exon skip
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
                        if (junction1[gene] == "e" and junction2[gene] == "s" and exon2[gene] - exon1[gene] > 1) or \
                           (junction2[gene] == "e" and junction1[gene] == "s" and exon1[gene] - exon2[gene] > 1): 
                            passGene.append(gene)

                            sc_size = gene2coding_size[gene] if gene in gene2coding_size else 0
                            # sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "Exon skipping"
                if len(inframe_gene) > 0: in_frame = "In-frame"

            ##########
            # check for alternative-3'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "e" and gene not in junction2 and gene in exon2) or \
                           (gene in junction2 and junction2[gene] == "e" and gene not in junction1 and gene in exon1):
                            passGene.append(gene)

                            sc_size = gene2coding_size[gene] if gene in gene2coding_size else 0
                            # sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "Alternative 3'SS"
                if len(inframe_gene) > 0: in_frame = "In-frame"


            # check for alternative-5'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "s" and gene not in junction2 and gene in exon2) or \
                           (gene in junction2 and junction2[gene] == "s" and gene not in junction1 and gene in exon1):
                            passGene.append(gene)

                            sc_size = gene2coding_size[gene] if gene in gene2coding_size else 0
                            # sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "Alternative 5'SS"
                if len(inframe_gene) > 0: in_frame = "In-frame"

            ##########
            # check for intronic-alternative-3'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "e" and gene not in junction2 and gene not in exon2) or \
                           (gene in junction2 and junction2[gene] == "e" and gene not in junction1 and gene not in exon1):
                            passGene.append(gene)
            
                if len(passGene) > 0: spliceClass = "Intronic alternative 3'SS"

            # check for intronic-alternative-5'-splice-site
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if (gene in junction1 and junction1[gene] == "s" and gene not in junction2 and gene not in exon2) or \
                           (gene in junction2 and junction2[gene] == "s" and gene not in junction1 and gene not in exon1):
                            passGene.append(gene)
                
                if len(passGene) > 0: spliceClass = "Intronic alternative 5'SS"

            ##########
            # within-exon
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2 and gene in exon1 and gene in exon2:
                        if exon1[gene] == exon2[gene]:
                            passGene.append(gene)

                            sc_size = gene2coding_size[gene] if gene in gene2coding_size else 0
                            # sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "Within exon"
                if len(inframe_gene) > 0: in_frame = "In-frame"

            ##########
            # check for exon-exon-junction
            if spliceClass == "":
                passGene = []
                inframe_gene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2:
                        if gene in exon1 and gene in exon2:
                            passGene.append(gene)

                            sc_size = gene2coding_size[gene] if gene in gene2coding_size else 0
                            # sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            if sc_size != 0 and sc_size % 3 == 0:
                                inframe_gene.append(gene)

                if len(passGene) > 0: spliceClass = "Exon exon junction"
                if len(inframe_gene) > 0: in_frame = "In-frame"

            ##########
            # check for within-gene 
            if spliceClass == "":
                passGene = []
                for gene in checkGenes:
                    if gene in gene1 and gene in gene2 and gene: passGene.append(gene)   
            
                if len(passGene) > 0: spliceClass = "Within gene"


            ##########
            # check for spliced-chimera 
            if spliceClass == "":
                passGene = []
                for g1 in gene1:
                    for g2 in gene2:
                        if (g1 in junction1 and junction1[g1] == "s" and g2 in junction2 and junction2[g2] == "e") or \
                           (g1 in junction1 and junction1[g1] == "e" and g2 in junction2 and junction2[g2] == "s"): 
                            passGene.append(g1 + ',' + g2)

                            # sc_size = gene2coding_size[gene] if gene in gene2coding_size else 0
                            # sc_size = spliced_coding_size(gene, None, chr_name, sj_start, sj_end, coding_tb, exon_margin)
                            # if sc_size != 0 and sc_size % 3 == 0:
                            #     inframe_gene.append(gene)


                if len(passGene) > 0: spliceClass = "Spliced chimera"
                # if len(inframe_gene) > 0: in_frame = "in-frame"

            ##########
            # check for unspliced-chimera 
            if spliceClass == "":
                passGene = []
                for g1 in gene1:
                    for g2 in gene2:
                        passGene.append(g1 + ',' + g2)

                if len(passGene) > 0: spliceClass = "Unspliced chimera"


            if spliceClass == "": spliceClass = "Other"
            

            # summarize the exon and junction information for display
            geneInfo1 = []
            exonInfo1 = []
            junctionInfo1 = []
            offsetInfo1 = []
            if len(gene1) > 0:
                for g1 in sorted(gene1):
                    if g1 not in passGene: continue 
                    geneInfo1.append(g1)
                    if g1 in exon1: 
                        exonInfo1.append(str(exon1[g1]))
                    else:
                        exonInfo1.append("*")

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
                exonInfo1.append("---")
                junctionInfo1.append("---")
                offsetInfo1.append("---")

            geneInfo2 = []
            exonInfo2 = []
            junctionInfo2 = []
            offsetInfo2 = []
            if len(gene2) > 0:
                for g2 in sorted(gene2):
                    if g2 not in passGene: continue
                    geneInfo2.append(g2)
                    if g2 in exon2:
                        exonInfo2.append(str(exon2[g2]))
                    else:
                        exonInfo2.append("*")
                    
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
                exonInfo2.append("---")
                junctionInfo2.append("---")
                offsetInfo2.append("---")

         

            print('\t'.join(F) + '\t' + spliceClass + '\t' + in_frame + '\t' + '\t'.join([';'.join(geneInfo1), ';'.join(exonInfo1), ';'.join(junctionInfo1), ';'.join(offsetInfo1), ';'.join(geneInfo2), ';'.join(exonInfo2), ';'.join(junctionInfo2), ';'.join(offsetInfo2)]), file = hout)
     
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc1.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc2.gene.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc1.exon.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc2.exon.bed"])
    subprocess.check_call(["rm", "-rf", output_file + ".tmp3.junc12.coding_size.bed"])

def juncmut_juncutils(input_file, output_file, control_list, genecode_gene_file, genome_id, is_grc, read_num_thres, junction_margin, exon_margin):
    import subprocess
    import shutil
    import os
    from junc_utils.utils import proc_star_junction

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

    juncutils_annotate(tmpfile2, output_file, genecode_gene_file, junction_margin, exon_margin, "gencode", genome_id, is_grc)

    for tfile in tmpfile_list:
        os.remove(tfile)

if __name__== "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    control_list = sys.argv[3].split(",")
    genecode_gene_file = sys.argv[4]

    juncmut_juncutils(input_file, output_file, control_list, genecode_gene_file, "hg38", "False", 1, 3, 30)
