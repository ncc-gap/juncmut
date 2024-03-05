#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_main(args):
    from .utils import check_reference
    from .juncmut_juncutils import juncmut_juncutils
    from .juncmut_assadj import juncmut_assadj
    from .juncmut_freq import juncmut_freq
    from .juncmut_mutpre import juncmut_mutpre
    from .juncmut_intersect import juncmut_intersect
    from .juncmut_rnamut import juncmut_rnamut
    from .juncmut_realign import juncmut_realign
    from .juncmut_supportread_count import juncmut_supportread_count
    from .juncmut_filt_bam import juncmut_filt_bam_main
    from .juncmut_annotgnomad import juncmut_annotgnomad

    genome_id, is_grc = check_reference(args.reference)

    juncmut_juncutils(
        args.input_file, 
        args.output_file+".fil.annot.txt", 
        args.control_file, 
        genome_id, 
        1
    )

    juncmut_assadj(
        args.output_file+".fil.annot.txt",
        args.output_file+".fil.annot.assadj.txt"
    )

    juncmut_freq(
        args.output_file+".fil.annot.assadj.txt", 
        args.output_file+".fil.annot.assadj.freq.txt",
        args.input_file, args.read_num_thres, args.freq_thres
    )

    juncmut_mutpre(
        args.output_file+".fil.annot.assadj.freq.txt",
        args.output_file+".fil.annot.assadj.freq.pmut.txt", 
        args.reference
    )

    juncmut_intersect(
        args.output_file+".fil.annot.assadj.freq.pmut.txt", 
        args.output_file+".fil.annot.assadj.freq.pmut.SJint.txt", 
        args.input_file
    )

    juncmut_rnamut(
        args.output_file+".fil.annot.assadj.freq.pmut.SJint.txt",
        args.output_file+".fil.annot.assadj.freq.pmut.SJint.rmut.txt", 
        args.rna_bam, 
        args.reference,
        args.mut_num_thres, 
        args.mut_freq_thres, 
    )

    juncmut_realign(
        args.output_file+".fil.annot.assadj.freq.pmut.SJint.rmut.txt",
        args.output_file+".fil.annot.assadj.freq.pmut.SJint.rmut.ed.txt", 
        args.rna_bam, 
        args.reference, 
        genome_id, 
        is_grc, 
        template_size = 10
    )

    juncmut_filt_bam_main(
        args.output_file+".fil.annot.assadj.freq.pmut.SJint.rmut.ed.txt",
        args.output_file, 
        args.rna_bam, 
        args.output_bam, 
        args.genecode_gene_file
    )

    if args.debug == "False":
        import os
        os.remove(args.output_file+".fil.annot.txt")
        os.remove(args.output_file+".fil.annot.assadj.txt")
        os.remove(args.output_file+".fil.annot.assadj.freq.txt")
        os.remove(args.output_file+".fil.annot.assadj.freq.pmut.txt")
        os.remove(args.output_file+".fil.annot.assadj.freq.pmut.SJint.txt")
        os.remove(args.output_file+".fil.annot.assadj.freq.pmut.SJint.rmut.txt")
        os.remove(args.output_file+".fil.annot.assadj.freq.pmut.SJint.rmut.ed.txt")

def annot_main(args):
    from .utils import check_reference
    from .juncmut_annotgnomad import juncmut_annotgnomad

    genome_id, is_grc = check_reference(args.reference)

    juncmut_annotgnomad(
        args.input_file,
        args.output_file,
        args.gnomad, 
        genome_id
    )

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
