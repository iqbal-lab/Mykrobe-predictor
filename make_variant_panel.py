#! /usr/bin/env python
from Bio import SeqIO
from mongoengine import connect
connect('atlas')
from vcf2db import VariantFreq
from panelgeneration import AlleleGenerator
from panelgeneration import Variant

kmer = 15
al = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta", kmer = kmer)


for vf in VariantFreq.objects().order_by('start'):
	context = [Variant(vft.reference_bases, vft.start , "/".join(vft.alternate_bases))  for vft in VariantFreq.objects(start__ne = vf.start, start__gt = vf.start - kmer, start__lt = vf.start + kmer)]
	variant = Variant(vf.reference_bases, vf.start , "/".join(vf.alternate_bases))
	panel = al.create(variant, context)
	for a in  panel.alts:
		print ">%s" % variant
		print a


