#! /usr/bin/env python

from .parser import create_parser


def main():

    cparser = create_parser()
    args = cparser.parse_args() #example 'program input_file --control_file --genome_id --reads_thres --freq_thres'
    if vars(args) == {}:
        cparser.print_usage()
    else:
        print(args)
        args.func(args)
            