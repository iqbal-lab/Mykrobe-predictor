#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import csv
import glob
## Read the kmer counts into a hash
import datetime

from atlas.vcf2db import CallSet
from atlas.vcf2db import GenotypedVariant
# from atlas.vcf2db import VariantPanel

from atlas.genotyping import ColourCovgsReader

import multiprocessing
import argparse
parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('coverage', metavar='coverage', type=str, help='File with coverage on alleles output from CORTEX align')
parser.add_argument('db_name', metavar='db_name', type=str, help='db_name')
parser.add_argument('--kmer', metavar='kmer', type=int, help='kmer size', default = 31)
parser.add_argument('--all', help='Store ref GT aswell as alt', default = False, action = "store_true")
args = parser.parse_args()

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

# Read fasta
try:
    call_set = CallSet.objects.get(name = args.sample)
except DoesNotExist:
    call_set = CallSet.create(name = args.sample)
## Clear any genotyped calls so far
GenotypedVariant.objects(call_set = call_set).delete()
gvs = []
with open(args.coverage, 'r') as infile:
    reader = ColourCovgsReader(infile)
    for allele in reader:
        allele_name, params = allele.name.split('?')
        alt_or_ref, name = allele_name.split('-')
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
        if alt_pnz:
            gvs.append(GenotypedVariant.create_object(name = name,
                                                      call_set = call_set,
                                                      ref_pnz = ref_pnz, 
                                                      alt_pnz = alt_pnz,
                                                      ref_coverage = ref_covg, 
                                                      alt_coverage = alt_covg,
                                                      gt = "1/1"))
        elif not alt_covg and args.all:
            gvs.append(GenotypedVariant.create_object(name = name,
                                                      call_set = call_set,
                                                      ref_pnz = ref_pnz, 
                                                      alt_pnz = alt_pnz,
                                                      ref_coverage = ref_covg, 
                                                      alt_coverage = alt_covg,
                                                      gt = "0/0")) 

# for vp in VariantPanel.objects():
#     alt_pnz, alt_covg = alt_coverage(vp.alts)
#     ref_pnz, ref_covg = coverage_on_seq(vp.ref)
#     print alt_pnz, alt_covg,  ref_pnz, ref_covg
#     if alt_pnz:
#         gvs.append(GenotypedVariant.create_object(name = vp.name,
#                                                   call_set = call_set,
#                                                   ref_pnz = ref_pnz, 
#                                                   alt_pnz = alt_pnz,
#                                                   ref_coverage = ref_covg, 
#                                                   alt_coverage = alt_covg,
#                                                   gt = "1/1"))
#     elif not alt_covg and args.all:
#         gvs.append(GenotypedVariant.create_object(name = vp.name,
#                                                   call_set = call_set,
#                                                   ref_pnz = ref_pnz, 
#                                                   alt_pnz = alt_pnz,
#                                                   ref_coverage = ref_covg, 
#                                                   alt_coverage = alt_covg,
#                                                   gt = "0/0"))        

GenotypedVariant.objects.insert(gvs)





