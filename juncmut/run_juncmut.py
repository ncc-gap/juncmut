#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:39:58 2019

@author: genome
"""

import os
from pathlib import Path

from .juncmut_juncutils import juncmut_juncutils 
from .juncmut_assadj import juncmut_assadj
from .juncmut_freq import juncmut_freq
from .juncmut_mutpre import juncmut_mutpre
from .juncmut_intersect import juncmut_intersect
from .juncmut_rnamut import juncmut_rnamut
from .juncmut_realign import juncmut_realign
from .juncmut_annotgnomad import juncmut_annotgnomad
from .utils import check_reference

def run_juncmut(args):
        
    genome_id, is_grc = check_reference(args.reference)
    pre = Path(args.input_file).stem.split('.')[0]
    os.makedirs("juncmut", exist_ok = True)

    juncmut_juncutils(args.input_file, "./juncmut/"+pre+".SJ.fil.annot.txt", args.control_file, genome_id, 1)

    juncmut_assadj("./juncmut/"+pre+".SJ.fil.annot.txt",
                   "./juncmut/"+pre+".SJ.fil.annot.assadj.txt")

    juncmut_freq("./juncmut/"+pre+".SJ.fil.annot.assadj.txt", 
                 "./juncmut/"+pre+".SJ.fil.annot.assadj.freq.txt",
                 args.input_file, args.read_num_thres, args.freq_thres)

    juncmut_mutpre("./juncmut/"+pre+".SJ.fil.annot.assadj.freq.txt",
                   "./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.txt", 
                   args.reference)

    juncmut_intersect("./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.txt", 
                     "./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.txt", args.input_file)

    juncmut_rnamut("./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.txt",
                   "./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.rmut.txt", 
                   args.rna_bam, args.reference)

    juncmut_realign("./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.rmut.txt",
                    "./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.rmut.ed.txt", 
                    args.rna_bam, args.reference, genome_id, is_grc, template_size = 10)

    juncmut_annotgnomad("./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.rmut.ed.txt",
                        "./juncmut/"+pre+".SJ.fil.annot.assadj.freq.pmut.SJint.rmut.ed.snp.txt",
                        args.gnomad, genome_id)

