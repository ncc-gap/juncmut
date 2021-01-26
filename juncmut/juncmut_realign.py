#! /usr/bin/env python3

"""
The program check bam or cram from the filename.
"""

import sys, random
from pathlib import Path
import pysam
import annot_utils
import edlib

def juncmut_realign(input_file, output_file, bam_file, reference, genome_id, is_grc, mut_num_thres, mut_freq_thres, template_size = 10):

    def generate_template_seq(template_file, ref_tb, junc_tb, mut_chr, mut_pos, mut_ref, mut_alt, 
                              junc_start, junc_end, junc_annotated, template_size):

        def change_base_check(tseq, tpos, tref, talt):
            # for debugging
            tseq_alt = tseq
            if tseq_alt[tpos] != tref:
                print("Inconsistent ref base", file = sys.stderr)
                sys.exit(1)
            tseq_alt = tseq[:tpos] + talt + tseq[(tpos + 1):]
            return(tseq_alt)

        hout = open(template_file, 'w') 

        ##########
        # no splicing
        # without mutation
        ir_ref_seq = ref_tb.fetch(region = mut_chr + ':' + str(mut_pos - template_size) + '-' + str(mut_pos + template_size - 1))
        print("No_SJ_Neg\t"+ir_ref_seq, file = hout)
        #print(ir_ref_seq, file = hout)

        # with mutation
        ir_alt_seq = change_base_check(ir_ref_seq, template_size, mut_ref, mut_alt)
        print("No_SJ_Pos\t"+ir_alt_seq, file = hout)
        #print(ir_alt_seq, file = hout)
        ##########

        ##########
        # target splicing
        # without mutation
        target_junc_ref_seq = ref_tb.fetch(region = mut_chr + ':' + str((junc_start - 1) - (template_size - 1)) + '-' + str(junc_start - 1)) + \
            ref_tb.fetch(region = mut_chr + ':' + str(junc_end + 1) + '-' + str(junc_end + 1 + (template_size - 1))) 
        print("Target_SJ_Neg\t"+target_junc_ref_seq, file = hout)
        #print(target_junc_ref_seq, file = hout)

        # with mutation
        if mut_pos >= (junc_start - 1) - (template_size - 1) and mut_pos <= junc_start - 1:
            target_junc_alt_seq = change_base_check(target_junc_ref_seq, mut_pos - ((junc_start - 1) - (template_size - 1)), mut_ref, mut_alt)
            print("Target_SJ_Pos\t"+target_junc_alt_seq, file = hout)
            #print(target_junc_alt_seq, file = hout)
        elif mut_pos >= junc_end + 1 and mut_pos <= (junc_end + 1) + (template_size - 1):
            target_junc_alt_seq = change_base_check(target_junc_ref_seq, mut_pos - ((junc_end + 1) + (template_size - 1)) - 1, mut_ref, mut_alt)
            print("Target_SJ_Pos\t"+target_junc_alt_seq, file = hout)
            #print(target_junc_alt_seq, file = hout)
        ##########

        ##########
        # normal splicing (multiple)
        
        records = junc_tb.fetch(region = mut_chr + ':' + str(mut_pos - template_size - 1) + '-' + str(mut_pos + template_size))
        njind = 0
        for record_line in records:
            record = record_line.split('\t')
            rj_start, rj_end = int(record[1]) - 1, int(record[2])
            if (junc_annotated == junc_start and rj_start == junc_start and abs(rj_end - mut_pos) <= template_size) or \
                (junc_annotated == junc_end and rj_end == junc_end and abs(rj_start - 1 - mut_pos) <= template_size):
                # without mutation
                normal_junc_ref_seq = ref_tb.fetch(region = mut_chr + ':' + str((rj_start - 1) - (template_size - 1)) + '-' + str(rj_start - 1)) + \
                    ref_tb.fetch(region = mut_chr + ':' + str(rj_end + 1) + '-' + str((rj_end + 1) + (template_size - 1)))
                print("Normal_SJ_Neg_" + str(njind) + "\t" + normal_junc_ref_seq, file = hout)
                #print(normal_junc_ref_seq, file = hout)
            
                if mut_pos >= (rj_start - 1) - (template_size - 1) and mut_pos <= rj_start - 1:
                    normal_junc_alt_seq = change_base_check(normal_junc_ref_seq, mut_pos - ((rj_start - 1) - (template_size - 1)), mut_ref, mut_alt)
                    print("Normal_SJ_Pos_" + str(njind) + "\t" + normal_junc_alt_seq, file = hout)  
                    #print(normal_junc_alt_seq)
                elif mut_pos >= rj_end + 1 and mut_pos <= (rj_end + 1) + (template_size - 1):
                    normal_junc_alt_seq = change_base_check(normal_junc_alt_seq, mut_pos - ((rj_end + 1) + (template_size - 1)) - 1, mut_ref, mut_alt)
                    print("Normal_SJ_Pos_" + str(njind) + "\t" + normal_junc_alt_seq, file = hout)
                    #print(normal_junc_alt_seq, file = hout)
            njind += 1
        ##########

        hout.close()


    def extract_read_around_boundary(bam_file, output_file, mut_chr, mut_pos, max_count = 10000):

        def check_read(read):

            check_flag = True 

            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip unmapped read  
            if flags[2] == "1" or flags[3] == "1": check_flag = False
     
            # skip supplementary alignmentx
            if flags[8] == "1" or flags[11] == "1": check_flag = False

            # skip duplicated reads
            if flags[10] == "1": check_flag = False

            return(check_flag)

        b_path = Path(bam_file)

        if b_path.suffix == '.bam':
            bamfile = pysam.AlignmentFile(bam_file, 'rb')
        if b_path.suffix == '.cram':
            bamfile = pysam.AlignmentFile(bam_file, 'rc')
          
        read_count = bamfile.count(region = mut_chr + ':' + str(mut_pos) + '-' + str(mut_pos), read_callback = check_read)
        read_ind_list = random.sample(range(read_count), min(max_count, read_count))
        read_inds = {}
        for i in read_ind_list:
            read_inds[i] = 1

        ind = -1
        hout = open(output_file, 'w') 
        for read in bamfile.fetch(region = mut_chr + ':' + str(mut_pos) + '-' + str(mut_pos)):
            
            if not check_read(read): continue

            ind = ind + 1
            if not ind in read_inds: continue

            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]
            #import pdb; pdb.set_trace()
            read_id = read.qname + '_1' if flags[6] == '1' else read.qname + '_2'
            read_id = read_id + '_' + str(read.reference_start + 1) + '_' + str(read.reference_end) + '_' + str(read.cigarstring) + '_' + str(read.next_reference_start + 1)
            print('>' + read_id + '\n' + read.seq, file = hout)

        bamfile.close()
        hout.close()
    """
    #remove duplicates
    def remove_duplicate_reads(input_fa, output_fa):
        read_dict = {}
        position_list =[]
        with open(input_fa, 'r') as fin:
            for row in fin:
                if row.startswith('>'):
                    read_id = row.rstrip('\n')
                    id_info = row.rstrip('\n').split('_')
                    position = id_info[1]+'_'+id_info[2]+'_'+id_info[4]
                    R = fin.readline().rstrip('\n')
                    if position in position_list: continue
                    else:
                        read_dict[read_id] = R
                        position_list.append(position)
                        
        with open(output_fa, 'w') as fout:
            for key, value in read_dict.items():
                fout.write(key + "\n")
                fout.write(value + "\n")
    """
    def process_cigar(res_cigar):
        #res_cigar = '16=1X1=1I1='
        # make a alignment object
        cigar_align = ''
        re_cigar = res_cigar.translate(str.maketrans({'=': ',=,','X': ',X,', 'I': ',I,', 'D': ',D,'}))
        cigar_list = re_cigar.split(',')
        del cigar_list[-1]
        for pnum in range(0,len(cigar_list),2):
            cigar_align = cigar_align + str(str(cigar_list[pnum+1]) * int(cigar_list[pnum]))
        proc_cigar_list = list(cigar_align)
        
        for i in range(0,5):
            pop_start = proc_cigar_list.pop(0) 
            pop_end = proc_cigar_list.pop(-1) 
            if pop_start == 'I' and pop_end == 'I':
                del proc_cigar_list[0]
                del proc_cigar_list[-1]
            elif pop_start == 'I' and pop_end != 'I':
                del proc_cigar_list[0]
            elif pop_start != 'I' and pop_end == 'I':
                del proc_cigar_list[-1]
            elif pop_start != 'I' and pop_end != 'I': continue
        
        mut_count = len([i for i in proc_cigar_list if i != '='])
        return(mut_count)    


    ref_tb = pysam.FastaFile(reference)
    #annot_utils.junction.make_junc_info(output_file + ".gencode.junc.bed.gz", "gencode", genome_id, is_grc, False)
    
    if str(is_grc) == 'False':
        annot_utils.junction.make_junc_info(output_file + ".gencode.junc.bed.gz", "gencode", genome_id, False, False)
    else:
        annot_utils.junction.make_junc_info(output_file + ".gencode.junc.bed.gz", "gencode", genome_id, True, False)
    junc_tb = pysam.TabixFile(output_file + ".gencode.junc.bed.gz")

    ## start 
    s_path = Path(input_file)
    sample_name = s_path.stem.split('.')[0]
    hout = open(output_file, 'w') 
    header = ["Mut_key", "SJ_key", "Sample", "SJ_Type", "SJ_Strand", "SJ_Read_Count", "SJ_Depth", "SJ_Freq",
              "Ref_Motif", "Possivle_Alt_Motif","Possible_Alt_key", "Is_GT/AG", "Is_in_exon","SJ_Overlap_Count", 
              "Chr","Mut_Pos", "Mut_Ref", "Mut_Alt", "Mut_Count", "Mut_Depth", "Mut_Freq",
              "Realign_No_SJ_Neg", "Realign_No_SJ_Pos", "Realign_Target_SJ_Neg", "Reaglin_Target_SJ_Pos",
              "Realign_Normal_SJ_Neg", "Realign_Normal_SJ_Pos","RNA_Mut"]
    print('\t'.join(header), file = hout)
    # for each mutation in SJ.
    with open(input_file, 'r') as hin:
        next(hin)
        for line in hin:
            # key_readid_dict to remove read duplicates
            F = line.rstrip('\n').split('\t')
            # Is a position of mutation in Exon or Intron
            #o-->
            if "5" in F[6] and "+" in F[7]: 
                dis_mut_SJstart = int(F[16])-int(F[1])
                if dis_mut_SJstart <0 and dis_mut_SJstart >-3:
                    is_exon = 'exon'
                elif dis_mut_SJstart <6 and dis_mut_SJstart >=0:
                    is_exon = 'intron'
                else:
                    is_exon = 'outside of motif'
            #-->o
            elif "3" in F[6] and "+" in F[7]: 
                dis_mut_SJend = int(F[16])-int(F[2])
                if dis_mut_SJend ==1:
                    is_exon = 'exon'
                elif dis_mut_SJend <1 and dis_mut_SJend >-6:
                    is_exon = 'intron'
                else:
                    is_exon = 'outside of motif'
            #o<--       
            elif "3" in F[6] and "-" in F[7]: 
                dis_mut_SJstart = int(F[16])-int(F[1])
                if dis_mut_SJstart ==-1:
                    is_exon = 'exon'
                elif dis_mut_SJstart <6 and dis_mut_SJstart >=0:
                    is_exon = 'intron'
                else:
                    is_exon = 'outside of motif'
            #<--o
            elif "5" in F[6] and "-" in F[7]: 
                dis_mut_SJend = int(F[16])-int(F[2])
                if dis_mut_SJend >0 and dis_mut_SJend <3:
                    is_exon = 'exon'
                elif dis_mut_SJend <1 and dis_mut_SJend >-6:
                    is_exon = 'intron'
                else:
                    is_exon = 'outside of motif'
            
            mut_key = F[0]+","+F[16]+","+F[17]+","+F[18]
            sj_key = F[0]+":"+F[1]+"-"+F[2]
            #if RNA_mutation reads>=2 and Freq>=0.05, do realign.
            if float(F[21]) < mut_freq_thres or int(F[20]) < mut_num_thres:
                print('\t'.join([mut_key, sj_key, sample_name, F[6], F[7], F[8], F[9], F[10], F[11], F[12], F[13], F[14],is_exon, F[15],
                                F[0], F[16], F[17], F[18], F[20], str(len(F[19])), F[21]]) +"\t-\t-\t-\t-\t-\t-\tFalse", file = hout)                
            else:
            # test if F[0]=='chr19' and str(F[16])=='11830133':  
                F = line.rstrip('\n').split('\t')
                mut_chr, junc_start, junc_end = F[0], int(F[1]), int(F[2])
                mut_pos, mut_ref, mut_alt = int(F[16]), F[17], F[18]
                if (F[6].endswith("5'SS") and F[7] == '+') or (F[6].endswith("3'SS") and F[7] == '-'):
                    junc_annotated = junc_end
                else:
                    junc_annotated = junc_start
    
                generate_template_seq(output_file + ".tmp.template.fa", ref_tb, junc_tb, mut_chr, mut_pos, mut_ref, mut_alt,
                    junc_start, junc_end, junc_annotated, template_size)
        
                extract_read_around_boundary(bam_file, output_file + ".tmp.read_seq.fa", mut_chr, mut_pos)
                
                #remove_duplicate_reads(output_file + ".tmp.read_seq.fa", output_file + ".tmp.read_seq_rmdup.fa")
    
                d_query = {}
                d_realign = {'No_SJ_Neg': 0, 'No_SJ_Pos': 0, 'Target_SJ_Neg': 0, 'Target_SJ_Pos': 0, 'Normal_SJ_Neg': 0, 'Normal_SJ_Pos': 0}            
                d_order = {'Normal_SJ_Neg': 0, 'No_SJ_Neg': 1, 'Target_SJ_Neg': 2, 'Normal_SJ_Pos': 3, 'No_SJ_Pos': 4, 'Target_SJ_Pos': 5}
                
                # query to dict    
                with open(output_file + ".tmp.template.fa", 'r') as qin:
                    for l in qin:
                        Q = l.rstrip('\n').split('\t')
                        d_query[Q[0]] = Q[1]
                
                # alignment
                with open(output_file + ".tmp.read_seq.fa", 'r') as rin:
                    rmut = '-'
                    
                    for row in rin:
                        
                        if row.startswith('>'):
                            #read_id = read_id + '_' + str(read.reference_start + 1) + '_' + str(read.reference_end) + '_' + str(read.cigarstring) + '_' + str(read.next_reference_start + 1) 
                            #read_id = ['>DRR016694.11533413', '2', '14065429', '14484178', '39M418647N15M2D47M', '14484206']

                            R = rin.readline().rstrip('\n')
                            
                            d_ed = {}
                            d_mut_count = {}
                            for key in d_query:
                                res = edlib.align(d_query[key],R, mode="HW", task="path")
                                res_cigar = res.get('cigar')
                                res_ed = res.get('editDistance')                             
                                #mut_count and edit distance
                                if res_ed <= 2:
                                    mut_count = process_cigar(res_cigar)
                                            
                                    # There may be multiple Normal SJ templates.  Choose template with smaller editDistance than 2.
                                    if key.startswith('Normal_SJ_Neg_'):
                                        if 'Normal_SJ_Neg' in d_ed:
                                            if d_ed['Normal_SJ_Neg'] <= res_ed: continue
                                            else:
                                                d_ed['Normal_SJ_Neg'] = res_ed
                                        else:
                                            d_ed['Normal_SJ_Neg'] = res_ed
                                        
                                        if 'Normal_SJ_Neg' in d_mut_count:
                                            if d_mut_count['Normal_SJ_Neg'] >= int(mut_count): continue
                                            else:
                                                d_mut_count['Normal_SJ_Neg'] = int(mut_count)
                                        else:
                                            d_mut_count['Normal_SJ_Neg'] = int(mut_count)
                                    
                                    elif key.startswith('Normal_SJ_Pos_'):
                                        if 'Normal_SJ_Pos' in d_ed:
                                            if d_ed['Normal_SJ_Pos'] <= res_ed: continue
                                            else:
                                                d_ed['Normal_SJ_Pos'] = res_ed
                                        else:
                                            d_ed['Normal_SJ_Pos'] = res_ed 

                                        if 'Normal_SJ_Pos' in d_mut_count:
                                            if d_mut_count['Normal_SJ_Pos'] >= int(mut_count): continue
                                            else:
                                                d_mut_count['Normal_SJ_Pos'] = int(mut_count)
                                        else:
                                            d_mut_count['Normal_SJ_Pos'] = int(mut_count)                                    
                    
                                    
                                    else:
                                        d_ed[key] = res_ed
                                        d_mut_count[key] = int(mut_count) 
                            
                            #choose minmum ed_distance and the key and count up.                                        
                            min_k = [kv[0] for kv in d_ed.items() if kv[1] == min(d_ed.values())]
                            if len(min_k) == 0: continue
                            else:
                                min_k_sort = sorted(min_k, key=lambda x: d_order[x])
                                #import pdb; pdb.set_trace()
                                # if multiple mutation for the selected key is exist, no count-up.
                                if int(d_mut_count[min_k_sort[0]]) == 0:
                                    new_count = d_realign[min_k_sort[0]]+1
                                    d_realign[min_k_sort[0]] = new_count  
                                    
                    if d_realign['Target_SJ_Pos'] + d_realign['No_SJ_Pos']+ d_realign['Normal_SJ_Pos'] >=1:
                        rmut = 'True'
                    else: rmut = 'False'                   
                    
                print('\t'.join([mut_key, sj_key, sample_name, F[6], F[7], F[8], F[9], F[10], F[11], F[12], F[13], F[14], is_exon, F[15],
                                         F[0], F[16], F[17], F[18], F[20], str(len(F[19])), F[21],
                                         str(d_realign["No_SJ_Neg"]), str(d_realign["No_SJ_Pos"]), 
                                         str(d_realign["Target_SJ_Neg"]), str(d_realign["Target_SJ_Pos"]),
                                         str(d_realign["Normal_SJ_Neg"]), str(d_realign["Normal_SJ_Pos"]), rmut]), file = hout)
                Path(output_file + ".tmp.template.fa").unlink()
                Path(output_file + ".tmp.read_seq.fa").unlink()
                

    ref_tb.close()
    junc_tb.close()
    
    Path(output_file + ".gencode.junc.bed.gz").unlink()
    Path(output_file + ".gencode.junc.bed.gz.tbi").unlink()

    hout.close()


if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("-input_file", metavar = "input_file", default = None, type = str,
                            help = "input file")         
    parser.add_argument("-output_file", metavar = "output_file", default = None, type = str,
                            help = "output files") 
    parser.add_argument("-bam_file", metavar = "bam_file", default = None, type = str,
                            help = "output files") 
    parser.add_argument("-reference", metavar = "reference", default = None, type = str,
                            help = "/full/path/to/reference")     
    parser.add_argument("-genome_id", metavar = "genome_id", default = None, type = str,
                            help = "genome_id")
    parser.add_argument("-is_grc", metavar = "is_grc", default = None, type = str,
                            help = "If chr prefix is in chr name, False") 
    parser.add_argument("-mut_num_thres", type = int, default = 2,
                        help = "If A mutation with the number of mutation alleles >= mut_num_thres and the frequency >= mut_freq_thres, do realign. (default: %(default)s)")    
    parser.add_argument("-mut_freq_thres", type = float, default = 0.05,
                        help = "If A mutation with the number of mutation alleles >= mut_num_thres and the frequency >= mut_freq_thres, do realign. (default: %(default)s)")
       
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    bam_file = args.bam_file
    reference = args.reference
    genome_id = args.genome_id
    is_grc = args.is_grc
    mut_num_thres = args.mut_num_thres 
    mut_freq_thres = args.mut_freq_thres
 
    juncmut_realign(input_file, output_file, bam_file, reference, genome_id, is_grc, mut_num_thres, mut_freq_thres, template_size = 10)
    


