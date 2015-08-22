#! /usr/bin/env python
import csv
from Bio import SeqIO
import glob
import numpy as np 
## Read the kmer counts into a hash
import sys
import datetime
from mongoengine import connect
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import ListField
from mongoengine import IntField
from mongoengine import ReferenceField
from vcf2db import CallSet
from vcf2db import GenotypedVariant
connect('atlas')

kmer_counts = {}
kmer = 15

for line in sys.stdin:
    line = line.rstrip('\n')
    row = line.split(' ')
    kmer_counts[row[0]] = int(row[1])




# Read fasta
def process_panel(filepath):
    coverage = {}
    for e,record in enumerate(SeqIO.parse(filepath,"fasta")):
        seq = str(record.seq)
        coverage_record = []
        for i in xrange(kmer+2):
            kmer_coverage = kmer_counts[seq[i:i+kmer]]
            if kmer_coverage <= 1:
                coverage_record = []
                break
            else:
                coverage_record.append(kmer_coverage)
        if coverage_record:
            coverage[record.name] = np.median(coverage_record)
    return coverage



def create_gvs(i):
    name, covg = i
    return 
call_set = CallSet.objects.get(name = "C00012831")
covg = process_panel("prelim_panel_k15.fasta")
gvs = []
for name, covg in covg.iteritems():
    gvs.append(GenotypedVariant.create_object(name = name, call_set = call_set, coverage = covg))
GenotypedVariant.objects.insert(gvs)


