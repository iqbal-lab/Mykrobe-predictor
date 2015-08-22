#! /usr/bin/env python
from Bio import SeqIO
from mongoengine import connect
connect('atlas')
from pymongo import MongoClient
client = MongoClient()
db = client.atlas
from vcf2db import VariantFreq
from vcf2db import Variant as CalledVariant
from vcf2db import VariantSet
from panelgeneration import AlleleGenerator
from panelgeneration import Variant

kmer = 15
al = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta", kmer = kmer)
# print "Extracting unique variants"
total_samples = VariantSet.objects.count()
results = CalledVariant.objects().aggregate({ "$group" : { "_id": { "name": "$name"}, "count": {"$sum": 1} }},
										 {"$match" : {"count" : {"$gt" : 0} } })
# print "Storing and sorting unique variants"
VariantFreq.drop_collection()
vfs = []
for result in results:
	var = CalledVariant.objects(name = result["_id"]["name"])[0]
	vf = VariantFreq.create(name = result["_id"]["name"],
					   count = result["count"],
					   total_samples = total_samples,
					   start = var['start'],
					   reference_bases = var['reference_bases'],
					   alternate_bases = var['alternate_bases']
					   )
	vfs.append(vf)
db.variant_freq.insert(vfs)



for vf in VariantFreq.objects().order_by('start'):
	context = [Variant(vft.reference_bases, vft.start , "/".join(vft.alternate_bases))  for vft in VariantFreq.objects(start__ne = vf.start, start__gt = vf.start - kmer, start__lt = vf.start + kmer)]
	variant = Variant(vf.reference_bases, vf.start , "/".join(vf.alternate_bases))
	panel = al.create(variant, context)
	for a in  panel.alts:
		print ">%s" % variant
		print a


