#! /usr/bin/env python
from __future__ import print_function


import os
import argparse
import sys
import atlas.version
import logging


def run_subtool(parser, args):
    if args.command == 'add':
        import atlas.commands.add as submodule
    elif args.command == "genotype":
        import atlas.commands.genotype as submodule
    elif args.command == "dump":
        import dump as submodule
    elif args.command == "amr":
        import atlas.commands.amr as submodule

    # run the chosen submodule.
    submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):

    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument(
            "-q",
            "--quiet",
            help="Do not output warnings to stderr",
            action="store_true",
            dest="quiet")

DEFAULT_KMER_SIZE = os.environ.get("KMER_SIZE", 21)
DEFAULT_DB_NAME = os.environ.get("DB_NAME", "atlas")


def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(
        prog='atlas',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", help="Installed atlas version",
                        action="version",
                        version="%(prog)s " + str(atlas.version.__version__))
    subparsers = parser.add_subparsers(
        title='[sub-commands]',
        dest='command',
        parser_class=ArgumentParserWithDefaults)
    #########################################
    # create the individual tool parsers
    #########################################

    ##########
    # Add
    ##########
    parser_add = subparsers.add_parser(
        'add',
        help='Adds a set of variants to the atlas')
    parser_add.add_argument('vcf', type=str, help='a vcf file')
    parser_add.add_argument('reference_set', type=str, help='reference set')
    parser_add.add_argument('-m','--method', type=str, help='variant caller method (e.g. CORTEX)', default = "NotSpecified")
    parser_add.add_argument('-f','--force',  action='store_true', help='Force recreate VariantSet')
    parser_add.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default=None)
    parser_add.set_defaults(func=run_subtool)

    # ##########
    # # Genotype
    # ##########
    parser_geno = subparsers.add_parser('genotype', help='Genotype a sample')
    parser_geno.add_argument(
        '-s',
        '--sample',
        type=str,
        help='sample id',
        required=True)
    parser_geno.add_argument(
        '-1',
        '--seq',
        type=str,
        help='Seq file',
        nargs='+',
        required=True)
    parser_geno.add_argument(
        '--panels',
        metavar='panels',
        type=str,
        nargs='+',
        help='panels',
        default=None)
    parser_geno.add_argument(
        '--name',
        metavar='name',
        type=str,
        help='name',
        default='atlas_gt')
    parser_geno.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default=None)
    parser_geno.add_argument(
        '-k',
        '--kmer',
        metavar='kmer',
        type=int,
        help='kmer size',
        default=None)
    parser_geno.add_argument(
        '--all',
        help='Store ref GT aswell as alt',
        default=False,
        action="store_true")
    parser_geno.add_argument(
        '--force',
        help='Force rebuilding of binaries',
        default=False,
        action="store_true")
    parser_geno.set_defaults(func=run_subtool)

    # ##########
    # # Dump panel
    # ##########
    parser_dump = subparsers.add_parser('dump',
                                        help='Dump a panel of variant alleles')
    parser_dump.add_argument(
        'ref',
        metavar='ref',
        type=str,
        help='reference filepath')
    parser_dump.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default=DEFAULT_DB_NAME)
    parser_dump.add_argument(
        '-k',
        '--kmer',
        metavar='kmer',
        type=int,
        help='kmer length',
        default=DEFAULT_KMER_SIZE)
    parser_dump.add_argument('--force', default=False, action="store_true")
    parser_dump.set_defaults(func=run_subtool)

    # ##########
    # # AMR predict
    # ##########
    parser_amr = subparsers.add_parser('amr',
                                       help="Predict the sample's antibiogram")
    parser_amr.add_argument(
        '-s',
        '--sample',
        type=str,
        help='sample id',
        required=True)
    parser_amr.add_argument(
        '-1',
        '--seq',
        type=str,
        help='Seq file',
        nargs='+',
        required=True)
    parser_amr.add_argument(
        '--panel',
        metavar='panel',
        type=str,
        help='panel',
        default=None)
    parser_amr.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default=None)
    parser_amr.add_argument(
        '-k',
        '--kmer',
        metavar='kmer',
        type=int,
        help='kmer length',
        default=DEFAULT_KMER_SIZE)
    parser_amr.add_argument(
        '--name',
        metavar='name',
        type=str,
        help='name',
        default='atlas_gt')
    parser_amr.add_argument(
        '--species',
        metavar='species',
        type=str,
        help='species',
        default=None)
    parser_amr.add_argument('--force', default=False, action="store_true")
    parser_amr.set_defaults(func=run_subtool)

    args = parser.parse_args()
    args.func(parser, args)


if __name__ == "__main__":
    main()
