#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:20:15 2019

@author: genome
"""

from .parser import create_parser


def main():

    cparser = create_parser()
    args = cparser.parse_args()
    if vars(args) == {}:
        cparser.print_usage()
    else:
        print(args)
        args.func(args)
