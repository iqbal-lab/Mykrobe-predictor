#! /usr/bin/env python
import datetime
import logging
logging.basicConfig(level=logging.DEBUG)
from Bio import SeqIO

from mongoengine import connect
connect('atlas')
from pymongo import MongoClient
import multiprocessing
client = MongoClient()
db = client.atlas
from vcf2db import VariantFreq
from vcf2db import Variant as CalledVariant
from vcf2db import VariantSet
from vcf2db import VariantPanel
from panelgeneration import AlleleGenerator
from panelgeneration import Variant
from utils import unique


kmer = 15
logging.info("start " + str(datetime.datetime.now() ) )

al = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta", kmer = kmer)
logging.info("Extracting unique variants " + str(datetime.datetime.now() ) )
## Get all variant freq names
current_vars = VariantFreq.objects().distinct('name')
new_names = CalledVariant.objects(name__nin = current_vars).distinct('name')
logging.info("%i new unique variants " % len(new_names) )



logging.info("Storing and sorting unique variants " + str(datetime.datetime.now() ) )
vfs = []
for name in new_names:
	reference_bases = name[0]
	alternate_bases = [name[-1]]
	start = int(name[1:-1])
	vf = VariantFreq.create(name = name,
					   # count = CalledVariant.objects(name = name).count(),
					   # total_samples = total_samples,
					   start = start,
					   reference_bases = reference_bases,
					   alternate_bases = alternate_bases
					   )
	vfs.append(vf)

logging.info("Inserting documents to DB " + str(datetime.datetime.now() ) )

def make_panel(vf):
	context = [Variant(vft.reference_bases, vft.start , "/".join(vft.alternate_bases)) for vft in VariantFreq.objects(start__ne = vf.start, start__gt = vf.start - kmer, start__lt = vf.start + kmer)]
	variant = Variant(vf.reference_bases, vf.start , "/".join(vf.alternate_bases))
	panel = al.create(variant, context)
	return VariantPanel().create_doc(vf, panel.alts)	

if vfs:
	db.variant_freq.insert(vfs)
	logging.info("Generating new panels " + str(datetime.datetime.now() ) )
	## Get names of panels that need updating 
	## Get all the variants that are within K bases of new variants
	update_names = []
	affected_variants = []
	for new_vf in VariantFreq.objects(name__in = new_names):
		query = VariantFreq.objects(start__gt = new_vf.start - kmer, start__lt = new_vf.start + kmer)
		for q in query:
			affected_variants.append(q.id)
			update_names.append(q.name)
	## Remove all panels that need updating
	VariantPanel.objects(variant__in = affected_variants).delete()
	## Make panels for all new variants and panels needing updating
	pool = multiprocessing.Pool(20)	
	variant_panels = pool.map(make_panel,  VariantFreq.objects(name__in = new_names + update_names ))
	new_panels = db.variant_panel.insert(variant_panels)

	logging.info("Writing panel " + str(datetime.datetime.now() ) )

	with open("panel_k15.fasta",'a') as panel_file:
		for variant_panel in VariantPanel.objects(id__in = new_panels):
			for a in variant_panel.alts:
				panel_file.write(">%s\n" % variant_panel.variant.name)
				panel_file.write("%s\n" % a)
	i = 0
	kmers = set()
	with open("panel_k15.kmers",'a') as panel_kmer_file:
		for a in VariantPanel.objects(id__in = new_panels).distinct('alts'):
			for i in xrange(len(a)-kmer+1):
				seq = a[i:i+kmer]
				if not seq in kmers:
					kmers.add(seq)
					panel_kmer_file.write(">%i\n" % i)
					panel_kmer_file.write("%s\n" % seq)				
					i+=1


logging.info("End " + str(datetime.datetime.now() ) )




