#! /usr/bin/env python 
import argparse
import sys
import atlas.version
import logging
logging.basicConfig(level=logging.ERROR)

def run_subtool(parser, args):
    if args.command == 'add':
        import atlas.commands.add as submodule
    elif args.command == "genotype":
        import atlas.commands.genotype as submodule
    elif args.command == "dump":
        import dump as submodule

    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

DEFAULT_DB_NAME = "atlas"
DEFAULT_KMER_SIZE = 31
def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='atlas', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", help="Installed atlas version",
                        action="version",
                        version="%(prog)s " + str(atlas.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)
    #########################################
    # create the individual tool parsers
    #########################################

    ##########
    # Add
    ##########
    parser_add = subparsers.add_parser('add',
                                        help='Adds a set of variants to the atlas')
    parser_add.add_argument('-s','--sample',  type=str, help='sample id')
    parser_add.add_argument('-f', '--vcf', metavar='vcf', type=str, help='a vcf file')
    parser_add.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default = None)
    parser_add.add_argument('--kmer', metavar='kmer', type=int, help='kmer length', default = None)
    parser_add.set_defaults(func=run_subtool)


    # ##########
    # # Genotype
    # ##########
    parser_geno = subparsers.add_parser('genotype', help='Genotype a sample')
    parser_geno.add_argument('-s','--sample',  type=str, help='sample id')
    parser_geno.add_argument('-1', '--seq', type=str, help='Seq file', nargs='+')
    parser_geno.add_argument('--name', metavar='name', type=str, help='name', default = 'atlas_gt')
    parser_geno.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default = None)
    parser_geno.add_argument('--kmer', metavar='kmer', type=int, help='kmer size', default = None)
    parser_geno.add_argument('--all', help='Store ref GT aswell as alt', default = False, action = "store_true")
    parser_geno.add_argument('--force', help='Force rebuilding of binaries', default = False, action = "store_true")
    parser_geno.set_defaults(func=run_subtool)

    # ##########
    # # Dump panel
    # ##########
    parser_dump = subparsers.add_parser('dump',
                                        help='Dump a panel of variant alleles')
    parser_dump.add_argument('ref', metavar='ref', type=str, help='reference filepath')    
    parser_dump.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default = DEFAULT_DB_NAME)
    parser_dump.add_argument('--kmer', metavar='kmer', type=int, help='kmer length', default = DEFAULT_KMER_SIZE)
    parser_dump.add_argument('--force', default = False, action = "store_true")
    parser_dump.set_defaults(func=run_subtool)

    args = parser.parse_args()
    args.func(parser, args)


if __name__ == "__main__":
    main()