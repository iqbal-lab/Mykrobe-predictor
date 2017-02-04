#! /usr/bin/env python
from __future__ import print_function

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
    for k, v in d.get("phylogenetics", {}).get("phylo_group", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_species_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("species", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_lineage_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("lineage", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_file_name(f):
    sample_name = os.path.basename(f).split('.')[0]
    return sample_name


def get_sample_name(f):
    return list(f.keys())[0]

def get_plate_name(f):
    return f.split('/')[-3]

def get_expected_depth(d):
    return str(d.get("expected_depth", -1))


def get_mean_read_length(d):
    return str(d.get("mean_read_length", -1))


def get_called_genes(d, drug=None):
    variants = []
    for variant_name, variant_call in d.items():
        if variant_call.get("_cls") == "Call.SequenceCall":
            per_cov = variant_call.get('info',{}).get('coverage',{}).get("percent_coverage")
            depth = variant_call.get('info',{}).get('coverage',{}).get("median_depth")
            variants.append(":".join([variant_name,
                             str(int(per_cov)),str(int(depth))]))
    return ";".join(variants)


def get_variant_calls(d):
    variants = []
    for variant_name, variant_call in d.items():
        if variant_call.get("_cls") != "Call.SequenceCall":
            wt_depth = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("median_depth")
            alt_depth = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("median_depth")

            wt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("percent_coverage")
            alt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("percent_coverage")
            if wt_per_cov < 100:
                wt_depth = 0
            if alt_per_cov <100:
                alt_depth =0 
            conf = variant_call.get('info',{}).get('conf',"1")

            variants.append(":".join([variant_name,
                             str(int(alt_depth)),str(int(wt_depth)), str(int(conf))
                           ]))
    return ";".join(variants)


if args.format == "long":
    header = [
        "mykrobe_version",
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
        "variants (gene:alt_depth:wt_depth:conf)",
        "genes (prot_mut-ref_mut:percent_covg:depth)"]
    print ("\t".join(header))
    rows = []
    for i, f in enumerate(args.files):
        file = get_file_name(f)
        try:
            d = load_json(f)
            k = list(d.keys())
            d = d[k[0]]
        except ValueError:
            d = {}
        

        phylo_group,phylo_group_per_covg,phylo_group_depth  = get_phylo_group_string(d)
        species,species_per_covg,species_depth  = get_species_string(d)
        lineage,lineage_per_covg,lineage_depth  = get_lineage_string(d)
        sample_name = get_sample_name(d)
        try:
            plate_name = get_plate_name(f)
        except  IndexError:
            plate_name = ""

        drug_list = sorted(d.get('susceptibility', {}).keys())
        drugs = sorted(drug_list)

        if not drugs:
            drugs = ["NA"]


        for drug in drugs:
            call = d.get('susceptibility', {}).get(drug, {})
            called_by_variants = get_variant_calls(call.get("called_by",{}))
            called_by_genes = get_called_genes(call.get("called_by",{}))
            row = [d.get("version",{}).get("mykrobe-predictor","-1"),
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
                call.get(
                    "predict",
                    'N'),
                called_by_variants, 
                called_by_genes
                ]
            # rows.append(row)
            print ("\t".join(row))

else:
    0 / 0
