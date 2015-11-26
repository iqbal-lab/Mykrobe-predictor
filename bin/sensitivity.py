#! /usr/bin/env python
## Compare genotype vs discovered
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import logging
logging.basicConfig(level=logging.DEBUG)
from atlas.vcf2db import *
import pymongo
import mongoengine
import argparse
from mongoengine import connect

from pymongo import MongoClient
client = MongoClient()




parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('db_name', metavar='db_name', type=str, help='db_name', default="atlas")
parser.add_argument('kmer', metavar='kmer', type=int, help='kmer size')
parser.add_argument('--diff', action = "store_true", help='show diff')
args = parser.parse_args()

DBNAME = 'atlas-%s-%i' % (args.db_name ,args.kmer) 
db = client[DBNAME]
connect(DBNAME)
# logging.info("Using DB %s"  % DBNAME)


## Get all discovered variants
variant_sets = VariantSet.objects(name__startswith = args.sample)
variants = Variant.objects(variant_set__in = variant_sets)#.order_by('start')
## All genotyped
discover_variant_set = set(variants.distinct('name_hash'))

atlas_gt_call_set_name = args.sample + "_atlas_gt"
atlas_gt_call_set = CallSet.objects.get(name = atlas_gt_call_set_name)
genotyped = GenotypedVariant.objects(call_set = atlas_gt_call_set, gt = "1/1")
genotyped_set = set(genotyped.distinct('name_hash'))

sample_variant_set = discover_variant_set.union(genotyped_set)

# print "CallSet", "union", "vars_in_callset", "intersection", "missing", "extra", "intersection/union"
print atlas_gt_call_set.name, len(sample_variant_set), len(genotyped_set), len(sample_variant_set & genotyped_set), len(sample_variant_set - genotyped_set), len(genotyped_set - sample_variant_set), float(len(sample_variant_set & genotyped_set))/float(len(sample_variant_set)) 
for variant_set in variant_sets:
    cur_variant_set = set(Variant.objects(variant_set = variant_set).distinct('name_hash'))
    print variant_set.name, len(sample_variant_set), len(cur_variant_set), len(sample_variant_set & cur_variant_set), len(sample_variant_set - cur_variant_set), len(cur_variant_set - sample_variant_set), float(len(sample_variant_set & cur_variant_set))/float(len(sample_variant_set)) 


# print "MISSING"
# for var in sample_variant_set - genotyped_set:
#     vp = VariantPanel.objects.get(name_hash = var)
#     v = vp.variant
#     print v.name

print GenotypedVariant.objects(name_hash__in = genotyped_set - discover_variant_set, call_set = atlas_gt_call_set, gt = "1/1").count()
if args.diff:
	print "EXTRA"
	for nh in genotyped_set - discover_variant_set:
	    # vp = VariantPanel.objects.get(name_hash = nh)
	    # v = vp.variant
	    # print v.name
	    gv = GenotypedVariant.objects.get(name_hash = nh, call_set = atlas_gt_call_set)

	    # if gv.gt == "1/1":
	    # 	print gv.to_json()
	    # 	print gv.alt_coverage

