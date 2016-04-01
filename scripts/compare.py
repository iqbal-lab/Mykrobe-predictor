#!/usr/bin/env python
# encoding: utf-8
"""
compare.py
Created on: 2016-03-30
Author: Phelim Bradley

Compares the JSON output of mykrobe predictor

"""
import sys

sys.path.append("/home/phelimb/git/Mykrobe-predictor")
import argparse
from mykrobe.utils import load_json
from mykrobe.utils import unique
from mykrobe.predict import MykrobePredictorSusceptibilityResult
parser = argparse.ArgumentParser(description='Compares the JSON output of mykrobe predictor')
parser.add_argument('truth', metavar='truth', type=str, help='truth')
parser.add_argument('--ana1', type=str, help='analyses', nargs='+', default = [])
parser.add_argument('--ana2', type=str, help='analyses', nargs='+', default = [])
args = parser.parse_args()

## First generate a table with a column for each commit for each sample for each drug
def file_paths_to_combined_dict(l):
	ana = {}
	for f in l:
		try:
			data = load_json(f)
		except ValueError, e:
			sys.stderr.write(str(e) + " %s \n" % f)
		else:
			assert data.keys()[0] not in ana
			ana.update(data)
	return ana

truth = load_json(args.truth)
ana1 = file_paths_to_combined_dict(args.ana1)
ana2 = file_paths_to_combined_dict(args.ana2)
def get_sample_ids(truth, ana1, ana2):
	return unique(truth.keys() + ana1.keys() + ana2.keys())

sample_ids = get_sample_ids(truth, ana1, ana2)
print len(sample_ids)

## sample  drug  truth  ana1   ana2 
## 1234    RIF   R     R     S 

## Report a summary of each analysis vs. truth

## ANA 1 
## Drug TP FP 
##

## Ana 2 

##
##


## Diff ana 1 ana 2 summary

##
##

## Report the diference between ana1 and and2

