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

def called_variant_list_to_panel_generated_variant_list(variants_in):
	variants_out = []
	for var in variants_in:
		for alt in var.alternate_bases:
			variants_out.append( Variant(var.reference_bases, var.start , alt) )
	return variants_out	

def seen_together(variants):
	## Takes a list of variants. 
	## Returns a list of variants that appear together (in the same variant set)
	## TODO - would be better if seen in same sample (accross several variant sets)
	variant_name_hashes = [v.name_hash for v in variants]
	variant_set_counter = Counter([v.variant_set.id for v in CalledVariant.objects(name_hash__in = variant_name_hashes)])
	variant_sets = [k for k,v in variant_set_counter.iteritems() if v > 1]
	contexts = []
	for vs in variant_sets:
		vars_together = [v for v in CalledVariant.objects(name_hash__in = variant_name_hashes, variant_set = vs)]
		if not vars_together in contexts:
			contexts.append(called_variant_list_to_panel_generated_variant_list(vars_together))
			variants = [var for var in variants if var not in vars_together] 
	for var in variants:
		contexts.append([var])
	return contexts + [[]]

def make_panels(vf):
	context = get_context(vf.start)
	panels = []
	contexts_seen_together = seen_together(context)	

	for alt in vf.alternate_bases:
		variant = Variant(vf.reference_bases, vf.start , alt)
		alts = []
		for context in contexts_seen_together:
			if len(context) <= 5:
				panel = al.create(variant, context)
				ref = panel.ref
				alts.extend(panel.alts)
		vp = VariantPanel().create(variant, vf, ref, alts)
		if not VariantPanel.objects(name_hash = vp.name_hash):
			panels.append(vp)
	return panels

def update_panel(vp):
	context = get_context(vp.variant.start)
	contexts_seen_together = seen_together(context)		
	
	ref, pos, alt = split_var_name(vp.name)
	variant = Variant(ref, pos, alt)

	alts = []
	for context in contexts_seen_together:
		if len(context) <= 5:
			panel = al.create(variant, context)
			alts.extend(panel.alts)
			ref = panel.ref
	vp = vp.update(ref, alts)
	return vp

def unique_name_hash(l):
	seen = set()
	return [seen.add(obj.name_hash) or obj for obj in l if obj.name_hash not in seen]	
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
	vfs = VariantFreq.objects.insert(vfs)	
	## Get names of panels that need updating 
	## Get all the variants that are within K bases of new variants
	## and that are not new variants
	print("Improving panel for genotyping" ) 
	update_name_hashes = []
	affected_variants = []
	if VariantPanel.objects().count() > 0:
		for new_vf in vfs:
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
	for vf in vfs:
		new_variant_panels.extend(make_panels(vf))
	for variant in affected_variants:
		vps= VariantPanel.objects(variant = variant)
		for vp in vps:
			updated_variant_panels.append(update_panel(vp))

	new_panels = VariantPanel.objects.insert(unique_name_hash(new_variant_panels)) + updated_variant_panels

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






