#! /usr/bin/env python

from __future__ import print_function
import gzip, subprocess
from . import chr_name, utils

def make_boundary_info(output_file, genome_id, is_grc, donor_size, acceptor_size):

    # create UCSC to GRC chr name corresponding table
    ucsc2grc = {}
    if is_grc:
        ucsc2grc = chr_name.make_ucsc2grc(genome_id)

    ucsc_gene_file = utils.set_ucsc_gene_file(genome_id, "refseq")


    donor_size_exon, donor_size_intron = [int(x) for x in donor_size.split(',')]
    acceptor_size_intron, acceptor_size_exon = [int(x) for x in acceptor_size.split(',')]

    key2junction, key2gene_id, key2exon_num = {}, {}, {}
    with gzip.open(ucsc_gene_file, 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if is_grc:
                if F[2] not in ucsc2grc: continue
                chr = ucsc2grc[F[2]]
            else:
                chr = F[2]
            # chr = ucsc2grc[F[2]] if F[2] in ucsc2grc else F[2] 

            starts = [int(x) for x in F[9].split(',') if x != '']
            ends = [int(x) for x in F[10].split(',') if x != '']
            strand = F[3]
            exon_num = int(F[8])
            gene_id = F[1]
            symbol = F[12]

            for i in range(0, exon_num - 1):
                if strand == '+': # donor
                    key = '\t'.join([chr, str(ends[i] - donor_size_exon), str(ends[i] + donor_size_intron), symbol, "donor", strand])
                else: # acceptor
                    key = '\t'.join([chr, str(ends[i] - acceptor_size_exon), str(ends[i] + acceptor_size_intron), symbol, "acceptor", strand])
                
                junction = chr + ':' + str(ends[i]) + '-' + str(starts[i + 1] + 1)

                if key not in key2junction: key2junction[key] = []
                if key not in key2gene_id: key2gene_id[key] = []
                if key not in key2exon_num: key2exon_num[key] = []

                key2junction[key].append(junction)
                key2gene_id[key].append(gene_id)
                key2exon_num[key].append(str(i))


            for i in range(1, exon_num):
                if strand == '+': # acceptor 
                    key = '\t'.join([chr, str(starts[i] - acceptor_size_intron), str(starts[i] + acceptor_size_exon), symbol, "acceptor", strand])
                else: # donor 
                    key = '\t'.join([chr, str(starts[i] - donor_size_intron), str(starts[i] + donor_size_exon), symbol, "donor", strand])

                junction = chr + ':' + str(ends[i - 1]) + '-' + str(int(starts[i]) + 1)

                if key not in key2junction: key2junction[key] = []
                if key not in key2gene_id: key2gene_id[key] = []
                if key not in key2exon_num: key2exon_num[key] = []
                
                key2junction[key].append(junction) 
                key2gene_id[key].append(gene_id) 
                key2exon_num[key].append(str(i)) 


    hout = open(output_file + ".unsorted.tmp", 'w')
    for key in sorted(key2junction):
        print('\t'.join([key, ','.join(key2junction[key]), ','.join(key2gene_id[key]), ','.join(key2exon_num[key])]), file = hout)
    hout.close()

    hout = open(output_file + ".sorted.tmp", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".unsorted.tmp"], stdout = hout)
    hout.close()

    hout = open(output_file, 'w')
    subprocess.check_call(["bgzip", "-f", "-c", output_file + ".sorted.tmp"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", output_file])

    subprocess.check_call(["rm", "-rf", output_file + ".unsorted.tmp"])
    subprocess.check_call(["rm", "-rf", output_file + ".sorted.tmp"])


