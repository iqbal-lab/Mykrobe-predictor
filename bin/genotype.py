#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import csv
import glob
import math
## Read the kmer counts into a hash
import datetime
from mongoengine import connect
from mongoengine import DoesNotExist
from atlas.vcf2db import CallSet
from atlas.vcf2db import GenotypedVariant
from atlas.vcf2db import VariantPanel

from atlas.genotyping import ColourCovgsReader

import multiprocessing
import argparse
parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('coverage', metavar='coverage', type=str, help='File with coverage on alleles output from CORTEX align')
parser.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default = 'tb')
parser.add_argument('--kmer', metavar='kmer', type=int, help='kmer size', default = 31)
parser.add_argument('--all', help='Store ref GT aswell as alt', default = False, action = "store_true")
args = parser.parse_args()

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

# Read fasta
try:
    call_set = CallSet.objects.get(name = args.sample + "_atlas_gt")
except DoesNotExist:
    call_set = CallSet.create(name = args.sample  + "_atlas_gt", sample_id = args.sample)
## Clear any genotyped calls so far
GenotypedVariant.objects(call_set = call_set).delete()
gvs = []

def max_pnz_threshold(vp):
    t =  max(100 - 2 * math.floor(float(max([len(alt) for alt in vp.alts])) / 100), 30)
    return t

with open(args.coverage, 'r') as infile:
    reader = ColourCovgsReader(infile)
    for allele in reader:
        allele_name, params = allele.name.split('?')
        alt_or_ref, _id = allele_name.split('-')
        vp = VariantPanel.objects.get(id = _id)
        MAX_PNZ_THRESHOLD = max_pnz_threshold(vp)
        if alt_or_ref == "ref":
            ref_pnz = allele.percent_non_zero_coverage
            ref_covg = allele.median_non_zero_coverage
            num_alts = int(params.split('=')[1])

            alt_pnz = 0
            alt_covg = 0            
            for _ in range(num_alts):
              allele = reader.next()
              if allele.percent_non_zero_coverage >= alt_pnz:
                  alt_pnz = allele.percent_non_zero_coverage
                  if allele.median_non_zero_coverage > alt_covg:
                      alt_covg = allele.median_non_zero_coverage 
        if alt_pnz >= MAX_PNZ_THRESHOLD and ref_pnz < MAX_PNZ_THRESHOLD:
            gt = "1/1"
        elif alt_pnz >= MAX_PNZ_THRESHOLD and ref_pnz >= MAX_PNZ_THRESHOLD:
            gt = "0/1"
        elif alt_covg < MAX_PNZ_THRESHOLD:
            gt = "0/0"
            
        if gt != "0/0" or args.all:
           gvs.append(GenotypedVariant.create_object(name = vp.name,
                                                      call_set = call_set,
                                                      ref_pnz = ref_pnz, 
                                                      alt_pnz = alt_pnz,
                                                      ref_coverage = ref_covg, 
                                                      alt_coverage = alt_covg,
                                                      gt = gt))

GenotypedVariant.objects.insert(gvs)





