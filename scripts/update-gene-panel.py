#! /usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
sys.path.append('/home/phelimb/git/atlas-core')
from atlas.utils import get_params

import json

import argparse
parser = argparse.ArgumentParser(description='Add length argument to panel')
parser.add_argument('dna_fasta', metavar='dna_fasta', type=str, help='dna_fasta')
parser.add_argument('json', metavar='json', type=str, nargs='+', help='json')
args = parser.parse_args()

seq_panel = {}
with open(args.dna_fasta, 'r') as infile:
	for i, record in enumerate(SeqIO.parse(infile, "fasta")):
		params = get_params(record.id)
		gene_name = params.get("name", i)
		version = int(params.get("version", i))
		if not gene_name in seq_panel:
			seq_panel[gene_name] = {}
		seq_panel[gene_name][version] = record.seq

for json_assem in args.json:
	with open(json_assem, 'r') as infile:
		try:
			data = json.load(infile)
		except ValueError:
			pass
		else:
			for gene_name, d in data.items():
				seq = d.get("dna").rstrip("*")
				if seq:
					if seq in seq_panel[gene_name].values():
						print("Already have version")
					else:
						print("Appending to panel")
						version = max(seq_panel[gene_name].keys()) + 1
						seq_panel[gene_name][version] = Seq(seq)

records = []
for gene_name, versions in seq_panel.items():
	for version, seq in versions.items():
		records.append(SeqRecord(seq,
			id = "%s?name=%s&version=%i" % (gene_name, gene_name, version) , description =""))
with open(args.dna_fasta, 'w') as outfile:
	SeqIO.write(records, outfile , 'fasta')




