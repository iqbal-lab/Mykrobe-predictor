#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import datetime
import logging
logging.basicConfig(level=logging.DEBUG)
from Bio import SeqIO

from mongoengine import connect
from pymongo import MongoClient
import multiprocessing
client = MongoClient()

from atlas.vcf2db import VariantFreq
from atlas.vcf2db import Variant as CalledVariant
from atlas.vcf2db import VariantSet
from atlas.vcf2db import VariantPanel
from atlas.vcf2db import split_var_name
from atlas.panelgeneration import AlleleGenerator
from atlas.panelgeneration import Variant

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

import argparse
parser = argparse.ArgumentParser(description='Parse VCF and upload variants to DB')
parser.add_argument('reference_filepath', metavar='reference_filepath', type=str, help='reference_filepath')
parser.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default="tb")
parser.add_argument('--kmer', metavar='kmer', type=int, help='kmer length', default = 31)
args = parser.parse_args()

db = client['atlas-%s-%i' % (args.db_name ,args.kmer) ]
connect('atlas-%s-%i' % (args.db_name ,args.kmer))

def make_panel(vf):
	context = [Variant(vft.reference_bases, vft.start , "/".join(vft.alternate_bases)) for vft in VariantFreq.objects(start__ne = vf.start, start__gt = vf.start - args.kmer, start__lt = vf.start + args.kmer)]
	variant = Variant(vf.reference_bases, vf.start , vf.alternate_bases)
	if len(context) <= 8:
		panel = al.create(variant, context)
		return VariantPanel().create_doc(vf, panel.ref, panel.alts)	

# "/home/phelimb/git/atlas/data/R00000022.fasta"
al = AlleleGenerator(reference_filepath = args.reference_filepath, kmer = args.kmer)
print("Extracting unique variants " )
## Get all variant freq names
current_vars = VariantFreq.objects().distinct('name_hash')
new_name_hashes = CalledVariant.objects(name_hash__nin = current_vars).distinct('name_hash')
print("%i new unique variants " % len(new_name_hashes) )
print("Storing and sorting unique variants " ) 
vfs = []
for name_hash in new_name_hashes:
	v = CalledVariant.objects(name_hash = name_hash)[0]
	vf = VariantFreq.create(name = v.name,
						name_hash = v.name_hash,
					   # count = CalledVariant.objects(name = name).count(),
					   # total_samples = total_samples,
					   start = v.start,
					   reference_bases = v.reference_bases,
					   alternate_bases = v.alternate_bases
					   )
	vfs.append(vf)

print("Inserting documents to DB " ) 

if vfs:
	VariantFreq.objects.insert(vfs)
	## Get names of panels that need updating 
	## Get all the variants that are within K bases of new variants
	print("Improving panel for genotyping" ) 
	update_name_hashes = []
	affected_variants = []
	if VariantPanel.objects().count() > 0:
		for new_vf in VariantFreq.objects(name_hash__in = new_name_hashes):
			query = VariantFreq.objects(start__gt = new_vf.start - args.kmer, start__lt = new_vf.start + args.kmer)
			for q in query:
				affected_variants.append(q.id)
				update_name_hashes.append(q.name)
		## Remove all panels that need updating
		VariantPanel.objects(variant__in = affected_variants).delete()
	## Make panels for all new variants and panels needing updating
	variant_panels = []
	for vf in VariantFreq.objects(name_hash__in = unique(new_name_hashes + update_name_hashes)):
		variant_panels.append(make_panel(vf))
	new_panels = db.variant_panel.insert(variant_panels)

	with open("panel_%s_k%i.fasta" % (args.db_name, args.kmer),'a') as panel_file:
		for variant_panel in VariantPanel.objects(id__in = new_panels):
			panel_file.write(">ref-%s?num_alts=%i\n" % (variant_panel.variant.id, len(variant_panel.alts)))
			panel_file.write("%s\n" % variant_panel.ref)
			for a in variant_panel.alts:
				panel_file.write(">alt-%s\n" % variant_panel.variant.id)
				panel_file.write("%s\n" % a)






