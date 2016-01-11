#! /usr/bin/env python 
from Bio import SeqIO


import argparse

parser = argparse.ArgumentParser(description='Add length argument to panel')
parser.add_argument('fasta_file', metavar='fasta_file', type=str, help='fasta_file')
args = parser.parse_args()

sequences = []
with open(args.fasta_file, 'r') as infile:
	for record in SeqIO.parse(infile, "fasta") :
		assert "?" in record.id
		assert "length" not in record.id
		record.id = record.id + "&length=%i" % len(record)
		sequences.append(record)


with open(args.fasta_file, 'w') as outfile:
	SeqIO.write(sequences, outfile, "fasta")

print "Done."


