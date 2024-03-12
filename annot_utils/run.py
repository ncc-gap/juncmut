#! /usr/bin/env python

from . import gene, exon, coding, junction, boundary, simple_repeat 

def gene_main(args):

    gene.make_gene_info(args.output_path, args.gene_model, args.genome_id, args.grc, args.add_ref_id)


def exon_main(args):

    exon.make_exon_info(args.output_path, args.gene_model, args.genome_id, args.grc, args.add_ref_id)


def coding_main(args):

    coding.make_coding_info(args.output_path, args.gene_model, args.genome_id, args.grc, args.add_ref_id)


def junction_main(args):

    junction.make_junc_info(args.output_path, args.gene_model, args.genome_id, args.grc, args.add_ref_id)


def boundary_main(args):

    boundary.make_boundary_info(args.output_path, args.genome_id, args.grc, args.donor_size, args.acceptor_size)


def simple_repeat_main(args):

    simple_repeat.make_simple_repeat_info(args.output_path, args.genome_id, args.grc)


