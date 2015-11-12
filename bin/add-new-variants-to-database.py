#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import logging
logging.basicConfig(level=logging.DEBUG)
from mongoengine import connect
from mongoengine import NotUniqueError
from mongoengine import OperationError

from atlas.vcf2db import CallSet
from atlas.vcf2db  import Reference
from atlas.vcf2db  import Variant
from atlas.vcf2db  import VariantSet
from atlas.vcf2db  import Call
from pymongo import MongoClient
client = MongoClient()


import vcf
import os
import csv

import argparse
parser = argparse.ArgumentParser(description='Parse VCF and upload variants to DB')
parser.add_argument('vcf', metavar='vcf', type=str, help='a vcf file')
parser.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default="tb")
parser.add_argument('--kmer', metavar='kmer', type=int, help='kmer length', default = 31)
args = parser.parse_args()

db = client['atlas-%s-%i' % (args.db_name ,args.kmer) ]
connect('atlas-%s-%i' % (args.db_name ,args.kmer))

def is_record_valid(record):
	valid = True
	for sample in record.samples:
		if sample["GT"] is None:
			valid = False
        if sample["GT_CONF"] < 1:
            valid = False
	return valid

vcf_reader = vcf.Reader(open(args.vcf, 'r'))
assert len(vcf_reader.samples) == 1

try:
	callset = CallSet.create(name = vcf_reader.samples[0], sample_id = vcf_reader.samples[0])
except NotUniqueError:
	callset = CallSet.objects.get(name = vcf_reader.samples[0])

try:
	variant_set = VariantSet.create(name = vcf_reader.samples[0])
except NotUniqueError:
	variant_set = VariantSet.objects.get(name = vcf_reader.samples[0])

try:
	reference = Reference.create(name = "R00000022")
except NotUniqueError:
	reference = Reference.objects.get(name = "R00000022")

variants = []
calls = []
for record in vcf_reader:
	if not record.FILTER and is_record_valid(record):
		for sample in record.samples:
			try:
				v = Variant.create_object(variant_set = variant_set,
				                    start = record.POS,
				                    reference_bases = record.REF,
								 	alternate_bases = [str(a) for a in record.ALT], 
								 	reference = reference)
				variants.append(v)
				c = Call.create_object(variant = v, call_set = callset, genotype = sample['GT'], genotype_likelihood = sample['GT_CONF'])
				calls.append(c)
			except (OperationError, ValueError) as e:
				print e

try:
	variant_ids =  [v.id for v in Variant.objects.insert(variants)]
except NotUniqueError:
	logging.error("Variants must be unique within a given VCF")
else:
	logging.info("Uploaded %i variants to database" % len(variant_ids))
	for i,_id in enumerate(variant_ids):
		calls[i]["variant"] = _id
	db.call.insert(calls)	

