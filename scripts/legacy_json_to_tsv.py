#! /usr/bin/env python

# This script is intended to load a JSON dict containing resistotypes,
# a list of comids and a list of drugs of interest. It will return a column for each drug,
# Where 1 = R, 0 S, 0 unknown.
import argparse
import json
import csv
import os
parser = argparse.ArgumentParser(description='''load a JSON dict containing resistotypes,
 a list of comids and a list of drugs of interest. It will return a column for each drug,
 Where 1 = R, 0 S, 0 unknown. ''')
parser.add_argument('--format', type=str,
                    help='--format', default="long")
parser.add_argument('files', type=str, nargs='+',
                    help='files')
args = parser.parse_args()


def load_json(f):
    with open(f, 'r') as infile:
        d = json.load(infile)
    return d


def get_drugs(drug_list):
    drugs = []
    for f in args.files:
        try:
            d = load_json(f)
        except ValueError:
            d = {}
        for drug in drug_list:
            if drug not in drugs:
                drugs.append(drug)
    return drugs


def get_phylo_group_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("phylo_group", {}).iteritems():
        s.append(k)
        per_cov.append(str(v))
    return ";".join(s), ".", ";".join(per_cov)


def get_species_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("species", {}).iteritems():
        s.append(k)
        per_cov.append(str(v))
    return ";".join(s), ".", ";".join(per_cov)



def get_lineage_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("lineage", {}).iteritems():
        s.append(k)
        per_cov.append(str(v))
    return ";".join(s), ".", ";".join(per_cov)


def get_file_name(f):
    sample_name = os.path.basename(f).split('.')[0]
    return sample_name


def get_sample_name(f):
    return f.split('/')[-2]

def get_plate_name(f):
    return f.split('/')[-3]

def get_expected_depth(d):
    return str(d.get("expected_depth", -1))


def get_mean_read_length(d):
    return str(d.get("mean_read_length", -1))


def get_called_genes(d, drug=None):
    genes = []
    for gene, coverage in d.get("called_genes", {}).iteritems():
        if coverage.get("induced_resistance") == drug:
            genes.append(":".join([gene,
                                   str(coverage.get("per_cov")),
                                   str(coverage.get('median_cov'))]))
    return ";".join(genes)

def get_called_variants(d, drug=None):
    variants = []
    for variant, coverage in d.get("called_variants", {}).iteritems():
        if coverage.get("induced_resistance") == drug:
            variants.append(":".join([variant,
                                   str(coverage.get("R_median_cov")),
                                   str(coverage.get('S_median_cov'))]))
    return ";".join(variants)


if args.format == "long":
    header = [
        "file",
        "plate_name",
        "sample",
        "drug",
        "phylo_group",
        "species",
        "lineage",
        "phylo_group_per_covg",
        "species_per_covg",
        "lineage_per_covg", 
        "phylo_group_depth",
        "species_depth",
        "lineage_depth",                  
        "susceptibility",
        "genes (gene:percent_coverage:depth)",
        "variants (prot_mut-ref_mut:alt_depth:wt_depth)"]
    print "\t".join(header)
    rows = []
    for i, f in enumerate(args.files):
        file = get_file_name(f)
        try:
            d = load_json(f)
        except ValueError:
            d = {}
        

        phylo_group,phylo_group_per_covg,phylo_group_depth  = get_phylo_group_string(d)
        species,species_per_covg,species_depth  = get_species_string(d)
        lineage,lineage_per_covg,lineage_depth  = get_lineage_string(d)
        sample_name = get_sample_name(f)
        plate_name = get_plate_name(f)

        drug_list = sorted(d.get('susceptibility', {}).keys())
        drugs = sorted(drug_list)

        if not drugs:
            drugs = ["NA"]


        for drug in drugs:
            call = d.get('susceptibility', {}).get(drug, "NA")
            called_by = get_called_genes(d, drug)
            called_by_variants  = get_called_variants(d, drug)
            row = [
                file,
                plate_name,
                sample_name,
                drug,
                phylo_group,
                species,
                lineage,
                phylo_group_per_covg,
                species_per_covg,
                lineage_per_covg,                  
                phylo_group_depth,
                species_depth,
                lineage_depth,                
                call,
                called_by,
                called_by_variants]
            print "\t".join(row)

else:
    0 / 0
