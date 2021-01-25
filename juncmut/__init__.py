__version__='0.4.1'

from .parser import create_parser


def main():

    cparser = create_parser()
    args = cparser.parse_args()
    if vars(args) == {}:
        cparser.print_usage()
    else:
        print(args)
        args.func(args)
