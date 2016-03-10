#! /usr/bin/env python
from __future__ import print_function
import sys
import os


import os
import argparse
import sys
sys.path.append(
    os.path.realpath(
        os.path.join(
            os.path.dirname(__file__),
            "..")))
from atlas.version import __version__
import logging

from base import ArgumentParserWithDefaults
from base import DEFAULT_KMER_SIZE


def run_subtool(parser, args):
    if args.command == "genotype":
        from atlas.cmds.genotype import run
    elif args.command == "predict":
        from atlas.cmds.amr import run

    # run the chosen submodule.
    run(parser, args)


def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(
        prog='mykrobe',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", help="mykrobe version",
                        action="version",
                        version="%(prog)s " + str(__version__))
    subparsers = parser.add_subparsers(
        title='[sub-commands]',
        dest='command',
        parser_class=ArgumentParserWithDefaults)

    sequence_parser_mixin = argparse.ArgumentParser(add_help=False)
    sequence_parser_mixin.add_argument(
        'sample',
        type=str,
        help='sample id')
    sequence_parser_mixin.add_argument(
        'seq',
        type=str,
        help='sequence files (fastq or bam)',
        nargs='+')
    sequence_parser_mixin.add_argument(
        '-k',
        '--kmer',
        metavar='kmer',
        type=int,
        help='kmer length (default:21)',
        default=DEFAULT_KMER_SIZE)
    sequence_parser_mixin.add_argument(
        '--tmp',
        help='tmp directory (default: /tmp/)',
        default="/tmp/")
    sequence_parser_mixin.add_argument(
        '--skeleton_dir',
        help='directory for skeleton binaries',
        default="atlas/data/skeletons/")
    sequence_parser_mixin.add_argument(
        '--mccortex31_path',
        help='Path to mccortex31',
        default="mccortex31")

    #########################################
    # create the individual tool parsers
    #########################################

    # ##########
    # # AMR predict
    # ##########

    parser_amr = subparsers.add_parser(
        'predict',
        parents=[sequence_parser_mixin],
        help="predict the sample's antibiogram")
    parser_amr.add_argument(
        'species',
        metavar='species',
        choices=['staph', 'tb'],
        type=str,
        help='species')
    parser_amr.add_argument(
        '--panel',
        metavar='panel',
        type=str,
        help='variant panel (default:bradley-2015)',
        choices=['bradley-2015', 'walker-2015'],
        default='bradley-2015')

    parser_amr.add_argument('--force', default=False, action="store_true")
    parser_amr.set_defaults(func=run_subtool)

    # ##########
    # # Genotype
    # ##########
    parser_geno = subparsers.add_parser(
        'genotype',
        parents=[sequence_parser_mixin],
        help='genotype a sample using a probe set')
    parser_geno.add_argument(
        'panels',
        metavar='panels',
        type=str,
        nargs='+',
        help='panels')
    parser_geno.add_argument(
        '--expected_depth',
        metavar='expected depth',
        type=int,
        help='expected depth',
        default=100)
    parser_geno.add_argument(
        '-f',
        '--force',
        help='Force rebuilding of binaries',
        default=False,
        action="store_true")
    parser_geno.set_defaults(func=run_subtool)

    args = parser.parse_args()
    args.func(parser, args)


if __name__ == "__main__":
    main()
