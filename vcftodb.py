#! /usr/bin/env python
from mongoengine import connect
from mongoengine import NotUniqueError
from mongoengine import OperationError

from vcf2db import CallSet
from vcf2db import Reference
from vcf2db import Variant
from vcf2db import VariantSet
from vcf2db import Call
from pymongo import MongoClient
client = MongoClient()
db = client.atlas
connect('atlas')

import vcf
import os
import csv

import argparse
parser = argparse.ArgumentParser(description='Parse VCF and upload variants to DB')
parser.add_argument('vcf', metavar='vcf', type=str, help='a vcf file')
args = parser.parse_args()

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
	variant_set = VariantSet.create(name = os.path.basename(args.vcf))
except NotUniqueError:
	variant_set = VariantSet.objects.get(name = os.path.basename(args.vcf))

try:
	reference = Reference.create(name = "R00000022")
except NotUniqueError:
	reference = Reference.objects.get(name = "R00000022")

variants = []
calls = []
for record in vcf_reader:
	if not record.FILTER and record.is_snp and is_record_valid(record):
		for sample in record.samples:
			try:
				v = Variant.create_object(variant_set = variant_set, start = record.POS, reference_bases = record.REF,
								 	alternate_bases = [str(a) for a in record.ALT], 
								 	reference = reference)
				variants.append(v)
				c = Call.create_object(variant = v, call_set = callset, genotype = sample['GT'], genotype_likelihood = sample['GT_CONF'])
				calls.append(c)
			except OperationError:
				pass


variant_ids =  db.variant.insert(variants)
for i,_id in enumerate(variant_ids):
	calls[i]["variant"] = _id
db.call.insert(calls)	

# import time
# time.sleep(20)
# Variant.ensure_indexes()
# Call.objects.insert(calls)			
