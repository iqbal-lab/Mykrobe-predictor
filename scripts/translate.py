#! /usr/bin/env python
from __future__ import print_function

from Bio import SeqIO

import argparse
parser = argparse.ArgumentParser(description='Add length argument to panel')
parser.add_argument('fasta_file', metavar='fasta_file', type=str, help='fasta_file')
parser.add_argument('out_fasta_file', metavar='out_fasta_file', type=str, help='out_fasta_file')
args = parser.parse_args()


sequences = []
proteins = set()

with open(args.fasta_file, 'r') as infile:
	for i, record in enumerate(SeqIO.parse(infile, "fasta")):
		record.seq = record.seq.translate()
		if str(record.seq) not in proteins:
			sequences.append(record)
			proteins.add(str(record.seq))

with open(args.out_fasta_file, 'w') as outfile:
	SeqIO.write(sequences, outfile, "fasta")
