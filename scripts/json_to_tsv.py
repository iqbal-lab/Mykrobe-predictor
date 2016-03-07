#! /usr/bin/env python

## This script is intended to load a JSON dict containing resistotypes, 
## a list of comids and a list of drugs of interest. It will return a column for each drug, 
## Where 1 = R, 0 S, 0 unknown. 
import argparse
import json
import csv
import os
parser = argparse.ArgumentParser(description='''load a JSON dict containing resistotypes, 
 a list of comids and a list of drugs of interest. It will return a column for each drug, 
 Where 1 = R, 0 S, 0 unknown. ''')
parser.add_argument('--format', type=str,
                   help='--format', default = "long")
parser.add_argument('files', type=str, nargs='+',
                   help='files')
args = parser.parse_args()

def load_json(f):
   with open(f, 'r') as infile:
      d =  json.load(infile)
   return d

def get_drugs(drug_list):
   drugs = []
   for f in args.files:
      try:
         d = load_json(f)
      except ValueError:
         d = {}
      for drug in drug_list:
         if not drug in drugs:
            drugs.append(drug)
   return drugs

def get_phylo_group_string(d):
   s = ""
   for k,v in d.get("phylogenetics",{}).get("phylo_group",{}).iteritems():
      s += "%s;" % (k) 
   return s

def get_species_string(d):
   s = ""
   for k,v in d.get("phylogenetics",{}).get("species",{}).iteritems():
      s += "%s;" % (k) 
   return s

def get_lineage_string(d):
   s = ""
   vmax = 0
   for k,v in d.get("phylogenetics",{}).get("lineage",{}).iteritems():
      if v > vmax:
         s += "%s;" % (k) 
   return s   

def get_file_name(f):
   sample_name = os.path.basename(f).split('.')[0]      
   return sample_name

def get_sample_name(f):
   return f.split('/')[-2]
def get_expected_depth(d):
   return str(d.get("expected_depth",-1))
def get_mean_read_length(d):
   return str(d.get("mean_read_length",-1))

def get_called_genes(d, drug = None):
   genes = []
   for gene, coverage in d.get("called_genes",{}).iteritems():
      if coverage.get("induced_resistance") == drug:
         genes.append(":".join([gene, str(coverage.get("per_cov")), str(coverage.get('median_cov'))]))
   return ";".join(genes)

def get_called_variants(d, drug = None):
   variants = []
   for gene, coverage in d.get("called_variants",{}).iteritems():
      if coverage.get("induced_resistance") == drug:
         variants.append(":".join([gene,  str(coverage.get('S_median_cov')),  str(coverage.get('R_median_cov'))]))
   return ";".join(variants)

if args.format == "long":
   header = ["file", "sample", "drug", "phylo_group","species", "lineage", "susceptibility", "variants"]   
   print "\t".join(header)
   rows = []
   for i,f in enumerate(args.files):
      try:
         d = load_json(f)
      except ValueError:
         d = {}

      phylo_group = get_phylo_group_string(d)
      species = get_species_string(d)
      lineage = get_lineage_string(d)
      file = get_file_name(f)
      sample_name = get_sample_name(f)

      drug_list = d[file].get('susceptibility',{}).keys()
      drug_list.sort()
      drugs = sorted(get_drugs(drug_list))

      if not drugs:
         drugs = ["NA"]
      for drug in drugs:
         called_genes = get_called_genes(d, drug = drug)
         called_variants = get_called_variants(d, drug = drug)
         call = d[file].get('susceptibility',{}).get(drug, {})
         row = [file, sample_name , drug, phylo_group, species,lineage,  call.get("predict", 'N'),  call.get("called_by","")]
         # rows.append(row)
         print "\t".join(row)

else:
   0/0




   



