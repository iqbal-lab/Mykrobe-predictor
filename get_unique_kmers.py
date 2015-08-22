#! /usr/bin/env python
from Bio import SeqIO
from utils import unique
kmer = 15
kmers = []
for record in SeqIO.parse("prelim_panel_k15.fasta","fasta"):
	a = str(record.seq)
	kmers.extend([a[i:i+kmer] for i in xrange(len(a)-kmer+1)])

for i, seq in enumerate(unique(kmers)):
	print ">%i" %i
	print seq

