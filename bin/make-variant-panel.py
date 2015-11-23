#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import datetime
import logging
from collections import Counter
logging.basicConfig(level=logging.DEBUG)
from Bio import SeqIO

from mongoengine import connect
import multiprocessing


from atlas.vcf2db import VariantFreq
from atlas.vcf2db import Variant as CalledVariant
from atlas.vcf2db import VariantSet
from atlas.vcf2db import VariantPanel
from atlas.vcf2db.models.base import split_var_name
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
parser.add_argument('--force', default = False, action = "store_true")
args = parser.parse_args()

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

def get_context(pos):
	context = []
	for vft in VariantFreq.objects(start__ne = pos, start__gt = pos - args.kmer, start__lt = pos + args.kmer):
		for alt in vft.alternate_bases:
			context.append( Variant(vft.reference_bases, vft.start , alt) )
	return context	

def make_panels(vf):
	context = get_context(vf.start)
	panels = []
	for alt in vf.alternate_bases:
		variant = Variant(vf.reference_bases, vf.start , alt)
		if len(context) <= 8:
			# print variant
			# print context
			panel = al.create(variant, context)
			panels.append(VariantPanel().create(variant,vf, panel.ref, panel.alts))
	return panels

def update_panel(vp):
	context = get_context(vp.variant.start)
	ref, pos, alt = split_var_name(vp.name)
	variant = Variant(ref, pos, alt) 
	if len(context) <= 8:
		panel = al.create(variant, context)
		vp = vp.update(panel.ref, panel.alts)
	return vp

# "/home/phelimb/git/atlas/data/R00000022.fasta"
al = AlleleGenerator(reference_filepath = args.reference_filepath, kmer = args.kmer)
print("Extracting unique variants " )
## Get all variant freq names
if args.force:
	new_name_hashes = CalledVariant.objects().distinct('name_hash')
	VariantFreq.objects().delete()
	VariantPanel.objects().delete()
else:
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



if vfs:
	print("Inserting %i documents to DB" % len(vfs)) 	
	VariantFreq.objects.insert(vfs)
	## Get names of panels that need updating 
	## Get all the variants that are within K bases of new variants
	## and that are not new variants
	print("Improving panel for genotyping" ) 
	update_name_hashes = []
	affected_variants = []
	if VariantPanel.objects().count() > 0:
		for new_vf in VariantFreq.objects(name_hash__in = new_name_hashes):
			query = VariantFreq.objects(start__gt = new_vf.start - args.kmer,
			                              start__lt = new_vf.start + args.kmer)
			for q in query:
				if not q.name_hash in new_name_hashes:
					affected_variants.append(q)
					update_name_hashes.append(q.name_hash)
		## Remove all panels that need updating
		# VariantPanel.objects(variant__in = affected_variants).delete()
	## Make panels for all new variants and panels needing updating
	new_variant_panels = []
	updated_variant_panels = []
	for vf in VariantFreq.objects(name_hash__in = new_name_hashes):
		new_variant_panels.extend(make_panels(vf))
	for variant in affected_variants:
		vps= VariantPanel.objects(variant = variant)
		for vp in vps:
			updated_variant_panels.append(update_panel(vp))

	new_panels = VariantPanel.objects.insert(new_variant_panels) + updated_variant_panels

	print("Appending %i panels" % len(new_panels) ) 
	if args.force:
		write_mode = "w"
	else:
		write_mode = "a"
	with open("panel_%s_k%i.fasta" % (args.db_name, args.kmer), write_mode) as panel_file:
		for variant_panel in VariantPanel.objects(id__in = [vp.id for vp in new_panels]):
			panel_file.write(">ref-%s?num_alts=%i\n" % (variant_panel.id, len(variant_panel.alts)))
			panel_file.write("%s\n" % variant_panel.ref)
			for a in variant_panel.alts:
				panel_file.write(">alt-%s\n" % variant_panel.id)
				panel_file.write("%s\n" % a)






