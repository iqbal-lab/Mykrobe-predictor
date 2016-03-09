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
from base import DEFAULT_DB_NAME
DEFAULT_KMER_SIZE = os.environ.get("KMER_SIZE", 31)

def run_subtool(parser, args):
    if args.command == 'add':
        from atlas.cmds.add import run
    elif args.command == "dump-probes":
        from atlas.cmds.dump import run
    elif args.command == "make-probes":
        from atlas.cmds.makeprobes import run

def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(
        prog='atlas',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", help="atlas version",
                        action="version",
                        version="%(prog)s " + str(__version__))
    subparsers = parser.add_subparsers(
        title='[sub-commands]',
        dest='command',
        parser_class=ArgumentParserWithDefaults)        


    db_parser_mixin = argparse.ArgumentParser(add_help=False)
    db_parser_mixin.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default=None)
    ##########
    # Add
    ##########
    parser_add = subparsers.add_parser(
        'add',
        help='adds a set of variants to the atlas',
        parents=[db_parser_mixin])
    parser_add.add_argument('vcf', type=str, help='a vcf file')
    parser_add.add_argument('reference_set', type=str, help='reference set')
    parser_add.add_argument(
        '-m',
        '--method',
        type=str,
        help='variant caller method (e.g. CORTEX)',
        default="NotSpecified")
    parser_add.add_argument(
        '-f',
        '--force',
        action='store_true',
        help='Force recreate VariantSet')
    parser_add.set_defaults(func=run_subtool)


    # ##########
    # # Dump panel
    # ##########
    parser_dump = subparsers.add_parser(
        'dump-probes',
        help='dump a panel of variant alleles',
        parents=[db_parser_mixin])
    parser_dump.add_argument(
        'reference_filepath',
        metavar='reference_filepath',
        type=str,
        help='reference_filepath')
    parser_dump.add_argument(
        '--kmer',
        metavar='kmer',
        type=int,
        help='kmer length',
        default=31)
    parser_dump.add_argument('--force', default=False, action="store_true")
    parser_dump.add_argument(
        '-v',
        '--verbose',
        default=False,
        action="store_true")
    parser_dump.set_defaults(func=run_subtool)

    ##################
    ### Make Probes ##
    ##################

    parser_make_probes = subparsers.add_parser(
        'make-probes', help='make probes from a list of variants',
        parents=[db_parser_mixin])
    parser_make_probes.add_argument(
        'reference_filepath',
        metavar='reference_filepath',
        type=str,
        help='reference_filepath')
    parser_make_probes.add_argument(
        '-v',
        '--variant',
        type=str,
        action='append',
        help='Variant in DNA positions e.g. A1234T',
        default=[])
    parser_make_probes.add_argument(
        '-f',
        '--file',
        type=str,
        help='File containing variants as rows A1234T')
    parser_make_probes.add_argument(
        '-g',
        '--genbank',
        type=str,
        help='Genbank file containing genes as features')
    parser_make_probes.add_argument(
        '-k',
        '--kmer',
        type=int,
        help='kmer length',
        default=31)
    parser_make_probes.set_defaults(func=run_subtool)   

    args = parser.parse_args()
    args.func(parser, args)


if __name__ == "__main__":
    main()
