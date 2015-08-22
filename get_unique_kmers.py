#! /usr/bin/env python
from Bio import SeqIO
from utils import unique
import argparse


parser = argparse.ArgumentParser(description='Parse VCF and upload variants to DB')
parser.add_argument('panel', metavar='panel', type=str, help='a fasta file')
args = parser.parse_args()

kmer = 15
kmers = []
for record in SeqIO.parse(args.panel,"fasta"):
	a = str(record.seq)
	kmers.extend([a[i:i+kmer] for i in xrange(len(a)-kmer+1)])

for i, seq in enumerate(unique(kmers)):
	print ">%i" %i
	print seq

