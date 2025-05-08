#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '0.6.1'

import sys
from .parser import create_parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    if vars(args) == {}:
        parser.print_usage()
        sys.exit(1)
    
    args.func(args)
