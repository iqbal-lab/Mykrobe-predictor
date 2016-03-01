#! /usr/bin/env python
from __future__ import print_function
import sys
import os


import os
import argparse
import sys
import atlas.version
import logging
sys.path.append("/..")


def run_subtool(parser, args):
    if args.command == 'add':
        from atlas.commands.add import run
    elif args.command == "genotype":
        from atlas.commands.genotype import run
    elif args.command == "dump-probes":
        from atlas.commands.dump import run
    elif args.command == "make-probes":
        from atlas.commands.makeprobes import run
    elif args.command == "predict":
        from atlas.commands.amr import run

    # run the chosen submodule.
    run(parser, args)


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
        prog='mykrobe',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", help="Installed mykrobe version",
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
        'panels',
        metavar='panels',
        type=str,
        nargs='+',
        help='panels')
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
        'kmer',
        metavar='kmer',
        type=int,
        help='kmer size')
    parser_geno.add_argument(
        '--all',
        help='Store ref GT aswell as alt',
        default=False,
        action="store_true")
    parser_geno.add_argument(
        '-f',
        '--force',
        help='Force rebuilding of binaries',
        default=False,
        action="store_true")
    parser_geno.add_argument(
        '--tmp',
        help='directory for skeleton binaries',
        default="/tmp/")
    parser_geno.set_defaults(func=run_subtool)

    # ##########
    # # Dump panel
    # ##########
    parser_dump = subparsers.add_parser('dump-probes',
                                        help='Dump a panel of variant alleles')
    parser_dump.add_argument(
        'reference_filepath',
        metavar='reference_filepath',
        type=str,
        help='reference_filepath')
    parser_dump.add_argument(
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default="tb")
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
        'make-probes', description='Make probes from a list of variants')
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
        '--db_name',
        metavar='db_name',
        type=str,
        help='db_name',
        default="tb")
    parser_make_probes.add_argument(
        '-k',
        '--kmer',
        type=int,
        help='kmer length',
        default=31)
    parser_make_probes.set_defaults(func=run_subtool)

    # ##########
    # # AMR predict
    # ##########
    parser_amr = subparsers.add_parser('predict',
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
