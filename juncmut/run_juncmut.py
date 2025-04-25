#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
def detect_main(args):
    #from .utils import check_reference
    from .juncmut_juncutils import juncmut_juncutils
    from .juncmut_assadj import juncmut_assadj
    from .juncmut_freq import juncmut_freq
    from .juncmut_mutpre import juncmut_mutpre
    from .juncmut_intersect import juncmut_intersect
    from .juncmut_rnamut import juncmut_rnamut
    from .juncmut_realign import juncmut_realign

    #genome_id, is_grc = check_reference(args.reference)

    juncmut_juncutils(
        args.sj_file, 
        args.output_file+".juncutils.txt", 
        args.control_file, 
        args.genecode_gene_file, 
        1, 3, 30, args.debug
    )

    juncmut_assadj(
        args.output_file+".juncutils.txt",
        args.output_file+".juncutils.assadj.txt"
    )

    juncmut_freq(
        args.output_file+".juncutils.assadj.txt", 
        args.output_file+".juncutils.assadj.freq.txt",
        args.sj_file, 
        args.read_num_thres, args.freq_thres
    )

    juncmut_mutpre(
        args.output_file+".juncutils.assadj.freq.txt",
        args.output_file+".juncutils.assadj.freq.mutpre.txt", 
        args.reference
    )

    juncmut_intersect(
        args.output_file+".juncutils.assadj.freq.mutpre.txt", 
        args.output_file+".juncutils.assadj.freq.mutpre.intersect.txt", 
        args.sj_file
    )

    juncmut_rnamut(
        args.output_file+".juncutils.assadj.freq.mutpre.intersect.txt",
        args.output_file+".juncutils.assadj.freq.mutpre.intersect.rnamut.txt", 
        args.bam_file, 
        args.reference,
        args.mut_num_thres, 
        args.mut_freq_thres, 
        args.support_read_rmdup_thres, 
    )

    juncmut_realign(
        args.output_file+".juncutils.assadj.freq.mutpre.intersect.rnamut.txt",
        args.output_file, 
        args.bam_file, 
        args.reference, 
        args.genecode_gene_file, 
        template_size = 10
    )

    if not args.debug:
        os.remove(args.output_file+".juncutils.txt")
        os.remove(args.output_file+".juncutils.assadj.txt")
        os.remove(args.output_file+".juncutils.assadj.freq.txt")
        os.remove(args.output_file+".juncutils.assadj.freq.mutpre.txt")
        os.remove(args.output_file+".juncutils.assadj.freq.mutpre.intersect.txt")
        os.remove(args.output_file+".juncutils.assadj.freq.mutpre.intersect.rnamut.txt")

def filt_bam_main(args):
    from .juncmut_filt_bam import juncmut_filt_bam

    juncmut_filt_bam(
        args.input_file, 
        args.input_bam, 
        args.output_bam, 
        args.genecode_gene_file
    )

def filt_main(args):
    from .juncmut_filt import juncmut_filt

    juncmut_filt(args.input_file, args.output_file)

def annot_main(args):
    import shutil
    from .utils import check_reference
    from .annot_gnomad import annot_gnomad
    from .annot_cgc import annot_cgc
    from .annot_clinvar import annot_clinvar
    from .annot_gene_list import annot_gene_list
    from .annot_cgd import annot_cgd

    if args.gnomad == "":
        shutil.copy(args.input_file, args.output_file + ".tmp1")
    else:
        genome_id, is_grc = check_reference(args.reference)
        annot_gnomad(args.input_file, args.output_file + ".tmp1", args.gnomad, genome_id)

    if args.cgc_file == "":
        shutil.copy(args.output_file + ".tmp1", args.output_file + ".tmp2")
    else:
        annot_cgc(args.output_file + ".tmp1", args.output_file + ".tmp2", args.cgc_file)

    if args.clinvar_file == "":
        shutil.copy(args.output_file + ".tmp2", args.output_file + ".tmp3")
    else:
        annot_clinvar(args.output_file + ".tmp2", args.output_file + ".tmp3", args.clinvar_file)

    if args.acmg_file == "":
        shutil.copy(args.output_file + ".tmp3", args.output_file + ".tmp4")
    else:
        annot_gene_list(args.output_file + ".tmp3", args.output_file + ".tmp4", args.acmg_file, "Is_ACMG")

    if args.clinvar_star234_file == "":
        shutil.copy(args.output_file + ".tmp4", args.output_file + ".tmp5")
    else:
        annot_gene_list(args.output_file + ".tmp4", args.output_file + ".tmp5", args.clinvar_star234_file, "Is_ClinVar_Star234")

    if args.pancan_file == "":
        shutil.copy(args.output_file + ".tmp5", args.output_file + ".tmp6")
    else:
        annot_gene_list(args.output_file + ".tmp5", args.output_file + ".tmp6", args.pancan_file, "Is_PancanAtlas")

    if args.dosage_sensitivity_file == "":
        shutil.copy(args.output_file + ".tmp6", args.output_file + ".tmp7")
    else:
        annot_gene_list(args.output_file + ".tmp6", args.output_file + ".tmp7", args.dosage_sensitivity_file, "Is_dosage_sensitivity")

    if args.cgd_file == "":
        shutil.copy(args.output_file + ".tmp7", args.output_file)
    else:
        annot_cgd(args.output_file + ".tmp7", args.output_file, args.cgd_file)

    if not args.debug:
        os.remove(args.output_file+".tmp1")
        os.remove(args.output_file+".tmp2")
        os.remove(args.output_file+".tmp3")
        os.remove(args.output_file+".tmp4")
        os.remove(args.output_file+".tmp5")
        os.remove(args.output_file+".tmp6")
        os.remove(args.output_file+".tmp7")

def sjclass_main(args):
    from .sjclass_transcript import sjclass_transcript
    from .sjclass_classify import sjclass_classify
    from .sjclass_frame import sjclass_frame

    sjclass_transcript(args.input_file, args.output_file + ".transcript.txt", args.gencode)
    sjclass_classify(args.output_file + ".transcript.txt", args.output_file + ".classify.txt", args.bam_file, args.sj_file, args.depth_th)
    sjclass_frame(args.output_file + ".classify.txt", args.output_file, args.reference)
    
    if not args.debug:
        os.remove(args.output_file +".transcript.txt")
        os.remove(args.output_file +".classify.txt")

"""
def validate_main(args):
    from .utils import check_reference
    from .juncmut_gmut import juncmut_gmut
    
    genome_id, is_grc = check_reference(args.reference)
    
    juncmut_gmut(
        args.input_file,
        args.output_file, 
        args.dna_bam, 
        args.reference, 
        is_grc, 
        args.mut_num_thres, 
        args.mut_freq_thres
    )
"""
