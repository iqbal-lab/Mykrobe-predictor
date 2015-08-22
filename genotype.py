#! /usr/bin/env python
import csv
from Bio import SeqIO
from multiprocessing.dummy import Pool as ThreadPool 
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
from mongoengine import FloatField

connect('atlas')

kmer_counts = {}
kmer = 15

for line in sys.stdin:
    line = line.rstrip('\n')
    row = line.split(' ')
    kmer_counts[row[0]] = int(row[1])

class GenotypedVariant(Document):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['name']
                }                                                   
                ]
            }    
    name = StringField()
    coverage = IntField()
    created_at = DateTimeField(required = True, default=datetime.datetime.now)
    
    start = IntField()
    reference_bases = StringField()
    alternate_bases = StringField()
    sample = StringField()

    @classmethod
    def create_object(cls, name, sample, coverage):
    	reference_bases = name[0]
    	start = int(name[1:-1])
    	alternate_bases = name[-1]
    	return cls( name = name, 
					reference_bases = reference_bases, 
					start = start, 
					alternate_bases = alternate_bases,
					sample = sample,
					coverage = int(coverage))
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

pool = ThreadPool(4) # Sets the pool size to 4
covg = process_panel("prelim_panel_k15.fasta")
gvs = []
for name, covg in covg.iteritems():
	gvs.append(GenotypedVariant.create_object(name = name, sample = "C00012831", coverage = covg))
GenotypedVariant.objects.insert(gvs)


