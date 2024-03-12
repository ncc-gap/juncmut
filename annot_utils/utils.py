#! /usr/bin/env python

import sys, pkg_resources, gzip, logging

from .logger import get_logger
logger = get_logger()

def set_ucsc_gene_file(genome_id, gene_model):
 
    if genome_id == "hg19":
        if gene_model == "refseq":
            ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/hg19/refGene.txt.gz")
        elif gene_model == "gencode":
            # ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/hg19/wgEncodeGencodeBasicV19.txt.gz")
            ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/hg19/wgEncodeGencodeCompV31lift37.txt.gz")
    elif genome_id == "hg38":
        if gene_model == "refseq":
            ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/hg38/refGene.txt.gz")
        elif gene_model == "gencode":
            # ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/hg38/wgEncodeGencodeBasicV24.txt.gz")
            ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/hg38/wgEncodeGencodeCompV31.txt.gz")
    elif genome_id == "mm10":
        if gene_model == "refseq":
            ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/mm10/refGene.txt.gz")
        elif gene_model == "gencode":
            # ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/mm10/wgEncodeGencodeBasicVM9.txt.gz")
            ucsc_gene_file = pkg_resources.resource_filename("annot_utils", "data/mm10/wgEncodeGencodeCompVM22.txt.gz")
    else:
        logger.error("genome_id shoud be hg19, hg38 or mm10.")
        sys.exit(1)

    return(ucsc_gene_file)


def grc_check(input_file, columns = [0], sep = '\t'):

    ucsc_chr_list = ['chr' + str(x) for x in list(range(1, 23)) + ['X', 'Y']]
    grc_chr_list = [str(x) for x in list(range(1, 23)) + ['X', 'Y']]
 
    input_chr_list = {}
    if input_file.endswith(".bam"):
        import pysam
        bamfile = pysam.AlignmentFile(input_file, 'rb')
        for tchr in bamfile.references:
            input_chr_list[tchr] = 1

    else:
        if input_file.endswith(".gz"):
            hin = gzip.open(input_file, 'rt')
        else:
            hin = open(input_file, 'r')

        for line in hin:
            F = line.rstrip('\n').split(sep)
            for col in columns:
                if col >= len(F): continue
                if F[col] not in input_chr_list: input_chr_list[F[col]] = 1

        hin.close()


    ucsc_flag = False
    grc_flag = False
    input_chr_list_ucsc = {}
    input_chr_list_grc = {}
    for tmp_chr in sorted(input_chr_list):
        if tmp_chr in ucsc_chr_list: 
            ucsc_flag = True
            input_chr_list_ucsc[tmp_chr] = 1

        if tmp_chr in grc_chr_list: 
            grc_flag = True
            input_chr_list_grc[tmp_chr] = 1
        
    if ucsc_flag == True and grc_flag == True:
        logger.error("UCSC (chr-prefixed) and GRC (non chr-prefixed) chromosome names are mixed! in " + input_file)
        logger.error("UCSC type chromosomes: " + ' '.join(input_chr_list_ucsc))
        logger.error("GRC type chromosomes: " + ' '.join(input_chr_list_grc))
        sys.exit(1)  
    elif ucsc_flag == True and grc_flag == False:
        return False
    elif ucsc_flag == False and grc_flag == True:
        return True
    else:
        logging.warning("No keys with valid chromosome names in " + input_file)
        return True


if __name__ == "__main__":
    import sys
    grc_check(sys.argv[1], [int(sys.argv[2])])

