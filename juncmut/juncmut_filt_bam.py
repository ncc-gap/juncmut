#! /usr/bin/env python


import subprocess, gzip

def juncmut_filt_bam_main(input_file, output_file, input_bam, output_bam, genecode_gene_file):

    def define_longest_transcript(splice_type, strand, tchr, normal_pos, ref_file):
        gene2tx_info = {}
        if splice_type == "5'SS" and strand == "+":
            with gzip.open(ref_file, 'rt') as hin:
                for record in hin:
                    R = record.rstrip('\n').split('\t')
                    ref_chr = R[2]
                    exonStarts = R[9].split(',')
                    #exonEnds = R[10].split(',')
                    if normal_pos in exonStarts and tchr == ref_chr:
                        ref_tx_id = R[1]
                        ref_tx_start = R[4]
                        ref_tx_end = R[5]
                        ref_gene = R[12]
                        gene2tx_info[ref_tx_id] = tchr, ref_tx_start, ref_tx_end, ref_gene
        if splice_type == "5'SS" and strand == "-":
            with gzip.open(ref_file, 'rt') as hin:
                for record in hin:
                    R = record.rstrip('\n').split('\t')
                    ref_chr = R[2]
                    #exonStarts = R[9].split(',')
                    exonEnds = R[10].split(',')
                    if normal_pos in exonEnds and tchr == ref_chr:
                        ref_tx_id = R[1]
                        ref_tx_start = R[4]
                        ref_tx_end = R[5]
                        ref_gene = R[12]
                        gene2tx_info[ref_tx_id] = tchr, ref_tx_start, ref_tx_end, ref_gene
        if splice_type == "3'SS" and strand == "+":
            with gzip.open(ref_file, 'rt') as hin:
                for record in hin:
                    R = record.rstrip('\n').split('\t')
                    ref_chr = R[2]
                    #exonStarts = R[9].split(',')
                    exonEnds = R[10].split(',')
                    if normal_pos in exonEnds and tchr == ref_chr:
                        ref_tx_id = R[1]
                        ref_tx_start = R[4]
                        ref_tx_end = R[5]
                        ref_gene = R[12]
                        gene2tx_info[ref_tx_id] = tchr, ref_tx_start, ref_tx_end, ref_gene
        if splice_type == "3'SS" and strand == "-":
            with gzip.open(ref_file, 'rt') as hin:
                for record in hin:
                    R = record.rstrip('\n').split('\t')
                    ref_chr = R[2]
                    exonStarts = R[9].split(',')
                    #exonEnds = R[10].split(',')
                    if normal_pos in exonStarts and tchr == ref_chr:
                        ref_tx_id = R[1]
                        ref_tx_start = R[4]
                        ref_tx_end = R[5]
                        ref_gene = R[12]
                        gene2tx_info[ref_tx_id] = tchr, ref_tx_start, ref_tx_end, ref_gene

        gene2chr = {}
        gene2start = {}
        gene2end = {}

        for id in gene2tx_info:
            tx_chr, tx_start, tx_end, gene = gene2tx_info[id]

            if gene in gene2chr:
                
                if gene2chr[gene] != tx_chr:
                    print('Wrong.')
            else: gene2chr[gene] = tx_chr
                        
            if gene in gene2start:
                if int(gene2start[gene]) > int(tx_start):
                    gene2start[gene] = tx_start
            else: gene2start[gene] = tx_start
            
            if gene in gene2end:
                if int(gene2end[gene]) < int(tx_end):
                    gene2end[gene] = tx_end
            else: gene2end[gene] = tx_end
        
        region_list=[]
        genes_list=[]
        genes = '---'
        for gene in gene2chr:
             
            region = gene2chr[gene] + ':' + str(gene2start[gene]) + '-' + str(gene2end[gene])
            region_list.append(region)
            genes_list.append(gene)
            genes = ','.join(genes_list)
            
        return(genes, region_list) 
    
    hout = open(output_file, 'w')
    ex_region_list =[]
    # open a file and make a trnscript list for RNA_Mut True.
    with open(input_file) as fin:
        header = fin.readline().rstrip('\n')
        new_header = header + "\tGene"
        print(new_header, file=hout)
        
        for line in fin:
            lie = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
            # col F[-3] is RNA_Mut
            mut_key = F[0].split(',')
            tchr = mut_key[0]
            strand = F[4]
            sj_pos = F[1].split(':')[1].split('-')
            sj_start = sj_pos[0]
            sj_end = sj_pos[1]
            
            #o-->
            if "5'SS" in F[3] and strand == "+": 
                splice_type = "5'SS"
                normal_pos = sj_end    
            #<--o
            if "5'SS" in F[3] and strand == "-":
                splice_type = "5'SS"
                normal_pos = int(sj_start)-1      
            #-->o
            if "3'SS" in F[3] and strand == "+":
                splice_type = "3'SS"
                normal_pos = int(sj_start)-1
            #o<--
            if "3'SS" in F[3] and strand == "-":
                splice_type = "3'SS"
                normal_pos = sj_end

            genes, region_list = define_longest_transcript(splice_type, strand, tchr, str(normal_pos), genecode_gene_file) 

            print(lie + '\t' + genes, file=hout)

            ex_region_list.extend(region_list)

    hout.close() 
           
    if not ex_region_list: 
    # initialize the file
        hout = open(output_bam + ".tmp.unsorted.sam", 'w')
        hout.close()
        #import pdb; pdb.set_trace() 
        hout = open(output_bam + ".tmp.unsorted.sam", 'a')
        
        for region in sorted(list(set(ex_region_list))):
            subprocess.check_call(["samtools", "view", input_bam, region], stdout = hout, stderr = subprocess.DEVNULL)
        hout.close()
    
        hout = open(output_bam + ".tmp.unsorted.rmdup.sam", 'w')
        subprocess.check_call(["sort", "-u", output_bam + ".tmp.unsorted.sam"], stdout = hout)
        hout.close()
        
        hout = open(output_bam + ".tmp.unsorted2.sam", 'w')
        subprocess.check_call(["samtools", "view", "-H", input_bam], stdout = hout, stderr = subprocess.DEVNULL)
        hout.close()
    
        hout = open(output_bam + ".tmp.unsorted2.sam", 'a')
        subprocess.check_call(["cat", output_bam + ".tmp.unsorted.rmdup.sam"], stdout = hout)
        hout.close()
    
        hout = open(output_bam + ".tmp.unsorted2.bam", 'w')
        subprocess.check_call(["samtools", "view", "-hbS", output_bam + ".tmp.unsorted2.sam"], stdout = hout)
        hout.close()
    
        hout = open(output_bam, 'w')
        subprocess.check_call(["samtools", "sort", output_bam + ".tmp.unsorted2.bam"], stdout = hout)
        hout.close()
    
        subprocess.check_call(["samtools", "index", output_bam])
            
        subprocess.check_call(["rm", "-rf", output_bam + ".tmp.unsorted.sam"])
        subprocess.check_call(["rm", "-rf", output_bam + ".tmp.unsorted.rmdup.sam"])
        subprocess.check_call(["rm", "-rf", output_bam + ".tmp.unsorted2.sam"])
        subprocess.check_call(["rm", "-rf", output_bam + ".tmp.unsorted2.bam"])
    

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("-input_file", metavar = "input_file", default = None, type = str,
                            help = "input file") 
    parser.add_argument("-output_file", metavar = "output_file", default = None, type = str,
                            help = "input file") 
    parser.add_argument("-input_bam", metavar = "input_bam", default = None, type = str,
                            help = "input bam") 
    parser.add_argument("-output_bam", metavar = "output_bam", default = None, type = str,
                            help = "output bam") 
    parser.add_argument("-gencode_gene_file", metavar = "gencode_gene_file", default = None, type = str,
                            help = "gencode_gene_file")  
    
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    input_bam = args.input_bam
    output_bam = args.output_bam
    gencode_gene_file = args.gencode_gene_file

    juncmut_filt_bam_main(input_file, output_file, input_bam, output_bam, gencode_gene_file)
    
