from ccreader import ColourCovgsReader
import argparse
parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('db_name', metavar='db_name', type=str, help='db_name')
parser.add_argument('kmer', metavar='kmer', type=int, help='kmer size')
parser.add_argument('kmer_count', metavar='kmer_count', type=str, help='kmer count')
parser.add_argument('--jobs', metavar='jobs', type=int, help='jobs', default = 10)
parser.add_argument('--all', help='Store ref GT aswell as alt', default = False, action = "store_true")
args = parser.parse_args()

reader = ColourCovgsReader()
