#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:39:58 2019

@author: genome
"""

from pathlib import Path

from .juncmut_env import juncmut_env
from .juncmut_juncutils import juncmut_juncutils 
from .juncmut_assadj import juncmut_assadj
from .juncmut_freq import juncmut_freq
from .juncmut_mutpre import juncmut_mutpre
from .juncmut_intersect import juncmut_intersect
from .juncmut_annotgnomadsnp import juncmut_annotgnomadsnp
from .juncmut_annotrnamut import juncmut_annotrnamut
from .juncmut_realign import juncmut_realign
from .utils import check_reference

def run_juncmut(args):
    
    # from . import juncmut_env, juncmut_juncutils, juncmut_assadj, juncmut_freq, juncmut_mutpre, juncmut_intersect, juncmut_annotgnomadsnp, juncmut_annotrnamut
    
    genome_id, is_grc = check_reference(args.reference)
    Path(args.output_file).parent.mkdir(parents = True, exist_ok = True)

    juncmut_juncutils(args.input_SJ, args.output_file + ".tmp.fil.annot.txt", args.control_file, genome_id, args.read_num_thres)
   
    juncmut_assadj(args.output_file + ".tmp.fil.annot.txt", 
                   args.output_file + ".tmp.fil.annot.assadj.txt")

    juncmut_freq(args.output_file + ".tmp.fil.annot.assadj.txt", 
                 args.output_file + ".tmp.fil.annot.assadjunifreqT.txt",
                 args.input_SJ, 
                 args.read_num_thres, args.freq_thres)

    juncmut_mutpre(args.output_file + ".tmp.fil.annot.assadjunifreqT.txt",
                   args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.txt", 
                   args.reference)

    juncmut_intersect(args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.txt", 
                      args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.SJinSJ.txt",
                      args.input_SJ)

    juncmut_annotrnamut(args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.SJinSJ.txt",
                        args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.SJinSJ.rmut.txt", 
                        args.rna_bam, args.reference)

    juncmut_annotgnomadsnp(args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.SJinSJ.rmut.txt",
                           args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.SJinSJ.rmut.snp.txt",
                           args.gnomad_path, genome_id)

    juncmut_realign(args.output_file + ".tmp.fil.annot.assadjunifreqT.pmut.SJinSJ.rmut.snp.txt",
                    args.output_file, 
                    args.rna_bam, args.reference, genome_id, is_grc, template_size = 10)

