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
# "/home/phelimb/git/atlas/data/R00000022.fasta"
al = AlleleGenerator(reference_filepath = args.reference_filepath, kmer = args.kmer)
print("Extracting unique variants " )
## Get all variant freq names
current_vars = VariantFreq.objects().distinct('name')
# result = VariantFreq._get_collection().aggregate([ 
#     { "$group": { "_id": "$name"}  }
#     ])
# current_vars = []
# for r in result:
# 	current_vars.append(r["_id"])
new_names = CalledVariant.objects(name__nin = current_vars).distinct('name')
# result = CalledVariant._get_collection().aggregate([ 
# 	{"$match" : {"name" : {"$nin" : current_vars}}},
#     { "$group": { "_id": "$name"}  }
#     ])

# new_names = []
# for r in result:
# 	new_names.append(r["_id"])

print("%i new unique variants " % len(new_names) )



print("Storing and sorting unique variants " ) 
vfs = []
for name in new_names:
	reference_bases, start, alternate_bases = split_var_name(name)
	vf = VariantFreq.create(name = name,
					   # count = CalledVariant.objects(name = name).count(),
					   # total_samples = total_samples,
					   start = start,
					   reference_bases = reference_bases,
					   alternate_bases = alternate_bases
					   )
	vfs.append(vf)

print("Inserting documents to DB " ) 

def make_panel(vf):
	context = [Variant(vft.reference_bases, vft.start , "/".join(vft.alternate_bases)) for vft in VariantFreq.objects(start__ne = vf.start, start__gt = vf.start - args.kmer, start__lt = vf.start + args.kmer)]
	variant = Variant(vf.reference_bases, vf.start , vf.alternate_bases)
	if len(context) <= 8:
		# print variant, context
		panel = al.create(variant, context)
		return VariantPanel().create_doc(vf, panel.ref, panel.alts)	

if vfs:

	db.variant_freq.insert(vfs)
	## Get names of panels that need updating 
	## Get all the variants that are within K bases of new variants
	print("Improving panel for genotyping" ) 
	update_names = []
	affected_variants = []
	if VariantPanel.objects().count() > 0:
		for new_vf in VariantFreq.objects(name__in = new_names):
			query = VariantFreq.objects(start__gt = new_vf.start - args.kmer, start__lt = new_vf.start + args.kmer)
			for q in query:
				affected_variants.append(q.id)
				update_names.append(q.name)
		## Remove all panels that need updating
		VariantPanel.objects(variant__in = affected_variants).delete()
	## Make panels for all new variants and panels needing updating
	# pool = multiprocessing.Pool(20)	
	# variant_panels = pool.map(make_panel,  VariantFreq.objects(name__in = unique(new_names + update_names) ))
	variant_panels = []
	for vf in VariantFreq.objects(name__in = unique(new_names + update_names)):
		variant_panels.append(make_panel(vf))
	new_panels = db.variant_panel.insert(variant_panels)


	

	kmers = set()

	with open("panel_%s_k%i.fasta" % (args.db_name, args.kmer),'a') as panel_file:
		for variant_panel in VariantPanel.objects(id__in = new_panels):
			panel_file.write(">ref-%s\n" % variant_panel.variant.name)
			panel_file.write("%s\n" % variant_panel.ref)
			for a in variant_panel.alts:
				panel_file.write(">alt-%s\n" % variant_panel.variant.name)
				panel_file.write("%s\n" % a)






