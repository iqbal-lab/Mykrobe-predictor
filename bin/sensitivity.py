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

call_set = CallSet.objects.get(name = args.sample)
## Get all discovered variants
variant_set = VariantSet.objects.get(name__startswith = args.sample)
variants = Variant.objects(variant_set = variant_set)#.order_by('start')
## All genotyped
genotyped = GenotypedVariant.objects(call_set = call_set, gt = "1/1")#.order_by('start')

variant_set = set(variants.distinct('name'))
genotyped_set = set(genotyped.distinct('name'))

print args.sample, len(variant_set), len(genotyped_set), len(variant_set & genotyped_set), len(variant_set - genotyped_set), len(genotyped_set - variant_set), float(len(variant_set & genotyped_set))/float(len(variant_set)) 

for var in variant_set - genotyped_set:
    print var, Variant.objects.get(name = var, id__in = [v.id for v in variants]).call.genotype_likelihood, VariantFreq.objects(name = var).count()
