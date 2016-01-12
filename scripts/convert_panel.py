#! /usr/bin/env python 
from Bio import SeqIO


import argparse

parser = argparse.ArgumentParser(description='Add length argument to panel')
parser.add_argument('fasta_file', metavar='fasta_file', type=str, help='fasta_file')
parser.add_argument('panel_type', metavar='panel_type', type=str, help='panel_type')
parser.add_argument('--name', metavar='name', type=str, help='name', default=None)
args = parser.parse_args()

if args.name is None:
	args.name = args.fasta_file.split(".")[0].split("/")[-1]

sequences = []
with open(args.fasta_file, 'r') as infile:
	for i, record in enumerate(SeqIO.parse(infile, "fasta")):
		assert "?" not in record.id
		record.id = record.id + "?name=%s&version=%i&length=%i&panel_type=%s" % (args.name, i + 1, len(record), args.panel_type)
		sequences.append(record)

with open(args.fasta_file, 'w') as outfile:
	SeqIO.write(sequences, outfile, "fasta")

print "Done."


