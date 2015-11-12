#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import csv
from Bio import SeqIO
import glob
from numpy import median
## Read the kmer counts into a hash
import datetime
from mongoengine import connect
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import ListField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import DoesNotExist
from atlas.vcf2db import CallSet
from atlas.vcf2db import GenotypedVariant
from atlas.vcf2db import VariantPanel

import multiprocessing
import argparse
parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('db_name', metavar='db_name', type=str, help='db_name')
parser.add_argument('kmer', metavar='kmer', type=int, help='kmer size')
parser.add_argument('kmer_count', metavar='kmer_count', type=str, help='kmer count')
parser.add_argument('--jobs', metavar='jobs', type=int, help='jobs', default = 10)
parser.add_argument('--all', help='Store ref GT aswell as alt', default = False, action = "store_true")
args = parser.parse_args()

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

kmer_counts = {}
kmer = args.kmer
with open(args.kmer_count,'r') as infile:
    reader = csv.reader(infile, delimiter = " ")
    for row in reader:
        try:
            kmer_counts[row[0]] = int(row[1])
        except IndexError:
            pass

# Read fasta
def alt_coverage(alts):
    best_percent_non_zero = 0
    best_coverage = 0
    for seq in alts:
        percent_non_zero, coverage = coverage_on_seq(seq)
        if percent_non_zero >= best_percent_non_zero:
            best_percent_non_zero = best_percent_non_zero
            if coverage > best_coverage:
                best_coverage = coverage
    return best_percent_non_zero, best_coverage

def coverage_on_seq(seq):
    coverage_record = []
    coverage = None
    for i in xrange(kmer+2):
        kmer_coverage = kmer_counts.get(seq[i:i+kmer],0)
        coverage_record.append(kmer_coverage)
    non_zero = [c for c in coverage_record if c > 0]
    percent_non_zero = sum(non_zero) / len(coverage_record)
    if non_zero:
        coverage = median(non_zero)
    return percent_non_zero, coverage  

try:
    call_set = CallSet.objects.get(name = args.sample)
except DoesNotExist:
    call_set = CallSet.create(name = args.sample)
## Clear any genotyped calls so far
GenotypedVariant.objects(call_set = call_set).delete()
gvs = []
for vp in VariantPanel.objects():
    alt_pnz, alt_covg = alt_coverage(vp.alts)
    ref_pnz, ref_covg = coverage_on_seq(vp.ref)
    if alt_covg:
        gvs.append(GenotypedVariant.create_object(name = vp.name,
                                                  call_set = call_set,
                                                  ref_pnz = ref_pnz, 
                                                  alt_pnz = alt_pnz,
                                                  ref_coverage = ref_covg, 
                                                  alt_coverage = alt_covg,
                                                  gt = "1/1"))
    elif not alt_covg and args.all:
        gvs.append(GenotypedVariant.create_object(name = vp.name,
                                                  call_set = call_set,
                                                  ref_pnz = ref_pnz, 
                                                  alt_pnz = alt_pnz,
                                                  ref_coverage = ref_covg, 
                                                  alt_coverage = alt_covg,
                                                  gt = "0/1"))        

GenotypedVariant.objects.insert(gvs)





