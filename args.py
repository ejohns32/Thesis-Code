import argparse
import os
import sys

sys.path.append(os.getcwd())


class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_arg_parser(description):
    return argparse.ArgumentParser(prog=sys.argv[0],
                                   description=description,
                                   formatter_class=SmartFormatter
                                   )
description = 'Stats printer for clustering algorithm.'
def_min_cluster = 5


def get_stats_args(description=description):
    arg_parser = get_arg_parser(description)
    arg_parser.add_argument('-f', '--filename',
                            action='store',
                            metavar='filename',
                            help='purity JSON file',
                            required=True)
    arg_parser.add_argument('-k', '--min-size', type=int,
                            default=def_min_cluster, help='minimum cluster size')
    return arg_parser
