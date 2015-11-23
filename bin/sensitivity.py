#! /usr/bin/env python
## Compare genotype vs discovered
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
from atlas.vcf2db import *
import pymongo
import mongoengine
import argparse
from mongoengine import connect

from pymongo import MongoClient
client = MongoClient()




parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('db_name', metavar='db_name', type=str, help='db_name', default="default")
parser.add_argument('kmer', metavar='kmer', type=int, help='kmer size')
args = parser.parse_args()

db = client['atlas-%s-%i' % (args.db_name ,args.kmer) ]
connect('atlas-%s-%i' % (args.db_name ,args.kmer))



## Get all discovered variants
variant_sets = VariantSet.objects(name__startswith = args.sample)
variants = Variant.objects(variant_set__in = variant_sets)#.order_by('start')
## All genotyped

sample_variant_set = set(variants.distinct('name_hash'))


atlas_gt_call_set_name = args.sample + "_atlas_gt"
atlas_gt_call_set = CallSet.objects.get(name = atlas_gt_call_set_name)
genotyped = GenotypedVariant.objects(call_set = atlas_gt_call_set, gt__ne = "0/0")
genotyped_set = set(genotyped.distinct('name_hash'))

print "CallSet", "union", "vars_in_callset", "intersection", "missing", "extra", "intersection/union"
print atlas_gt_call_set.name, len(sample_variant_set), len(genotyped_set), len(sample_variant_set & genotyped_set), len(sample_variant_set - genotyped_set), len(genotyped_set - sample_variant_set), float(len(sample_variant_set & genotyped_set))/float(len(sample_variant_set)) 
for variant_set in variant_sets:
    cur_variant_set = set(Variant.objects(variant_set = variant_set).distinct('name_hash'))
    print variant_set.name, len(sample_variant_set), len(cur_variant_set), len(sample_variant_set & cur_variant_set), len(sample_variant_set - cur_variant_set), len(cur_variant_set - sample_variant_set), float(len(sample_variant_set & cur_variant_set))/float(len(sample_variant_set)) 


# print "MISSING"
# for var in sample_variant_set - genotyped_set:
#     vp = VariantPanel.objects.get(name_hash = var)
#     v = vp.variant
#     print v.name

# print "EXTRA"
# for var in genotyped_set - sample_variant_set:
#     print var
#     vp = VariantPanel.objects.get(name_hash = var)
#     v = vp.variant
#     print v.name

