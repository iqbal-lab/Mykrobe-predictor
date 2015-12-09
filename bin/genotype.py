#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import csv
import glob
import math
import subprocess
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
parser.add_argument('-s','--sample',  type=str, help='sample id')
parser.add_argument('-1', '--seq', type=str, help='Seq file', nargs='+')
parser.add_argument('--name', metavar='name', type=str, help='name', default = 'atlas_gt')
parser.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default = None)
parser.add_argument('--kmer', metavar='kmer', type=int, help='kmer size', default = None)
parser.add_argument('--all', help='Store ref GT aswell as alt', default = False, action = "store_true")
parser.add_argument('--force', help='Force rebuilding of binaries', default = False, action = "store_true")
args = parser.parse_args()

if args.db_name is None:
    args.db_name = os.environ.get("DB_NAME")
if args.db_name is None:
    raise ValueError("db_name needs to be set. Either run with --db_name :db_name or export DB_NAME=:db_name")

if args.kmer is None:
    args.kmer = os.environ.get("KMER_SIZE")
if args.kmer is None:
    raise ValueError("kmer needs to be set. Either run with --kmer :kmer_size or export KMER_SIZE=:kmer_size")
else:
    args.kmer = int(args.kmer)

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

## If panel does not exists then build it
panel_filepath = os.path.abspath("data/panels/panel-%s-%i.fasta" % (args.db_name, args.kmer)) 
if not os.path.exists(panel_filepath):
    raise ValueError("Could not find a panel at %s. Run 'atlas dump'. " % panel_filepath)
## If ctx binary does not exist then build it
ctx_skeleton_filepath = os.path.abspath("data/skeletons/panel-%s-%i.ctx" % (args.db_name, args.kmer)) 
if not os.path.exists(ctx_skeleton_filepath):
    subprocess.check_output(["/home/phelimb/git/mccortex/bin/mccortex31", "build", "-q",
                             "-k", str(args.kmer), "-s", "panel-%s-%i" % (args.db_name, args.kmer),
                             "-1", panel_filepath, ctx_skeleton_filepath])
## Now get coverage on panel
sample = "-".join([args.sample, args.db_name, str(args.kmer)])
ctx_tmp_filepath = "/tmp/%s.ctx" % sample
covg_tmp_file_path = "/tmp/%s.covg" % sample

if not os.path.exists(ctx_tmp_filepath) or not os.path.exists(covg_tmp_file_path) or args.force:
    if os.path.exists(ctx_tmp_filepath):
        os.remove(ctx_tmp_filepath)
    if os.path.exists(covg_tmp_file_path):
        os.remove(covg_tmp_file_path)      
    cmd = ["/home/phelimb/git/mccortex/bin/mccortex31", "geno", "-q",
           "-I", ctx_skeleton_filepath,
           "-k", str(args.kmer), "-s", sample,
           "-o", covg_tmp_file_path]
    for seq in args.seq:
      cmd.extend(["-1", seq])
    cmd.extend(["-c", panel_filepath, "/tmp/%s.ctx" % sample])
    subprocess.check_output(cmd)



try:
    call_set = CallSet.objects.get(name = args.sample + "_%s" % args.name)
except DoesNotExist:
    call_set = CallSet.create(name = args.sample  + "_%s" % args.name, sample_id = args.sample)
## Clear any genotyped calls so far
GenotypedVariant.objects(call_set = call_set).delete()
gvs = {}

def max_pnz_threshold(vp):
    t =  max(100 - 2 * math.floor(float(max([len(alt) for alt in vp.alts])) / 100), 30)
    return t

def get_params(url):
    params = {}
    p_str = url.split("?")[1]
    p_str = p_str.split('&')
    for p in p_str:
        k,v = p.split("=")
        params[k] = v
    return params

with open("/tmp/%s.covg" % sample, 'r') as infile:
    reader = csv.reader(infile, delimiter = "\t")
    for row in reader:
        allele = row[0]
        median_coverage = int(row[1])
        pnz = float(row[2])

        allele_name = allele.split('?')[0]
        params = get_params(allele)
        alt_or_ref, _id = allele_name.split('-')
        # vp = VariantPanel.objects.get(id = _id)
        MAX_PNZ_THRESHOLD = 1.0#max_pnz_threshold(vp)
        if alt_or_ref == "ref":
            ref_pnz = pnz
            ref_covg = median_coverage

            num_alts = int(params.get("num_alts"))
            alt_pnz = 0
            alt_covg = 0            
            for _ in range(num_alts):
              row = reader.next()
              allele = row[0]
              median_coverage = int(row[1])
              pnz = float(row[2])      
              if pnz >= alt_pnz:
                  alt_pnz = pnz
                  if median_coverage > alt_covg:
                      alt_covg = median_coverage          
        if alt_pnz >= MAX_PNZ_THRESHOLD and ref_pnz < MAX_PNZ_THRESHOLD and alt_covg > 0:
            gt = "1/1"
        elif alt_pnz >= MAX_PNZ_THRESHOLD and ref_pnz >= MAX_PNZ_THRESHOLD and alt_covg > 0:
            gt = "0/1"
        elif alt_covg < MAX_PNZ_THRESHOLD:
            gt = "0/0"
        if gt != "0/0" or args.all:
           print allele_name, params.get("gene"), params.get("mut")
           gvs[_id] = GenotypedVariant.create_object(name = _id,
                                                      call_set = call_set,
                                                      ref_pnz = ref_pnz, 
                                                      alt_pnz = alt_pnz,
                                                      ref_coverage = ref_covg, 
                                                      alt_coverage = alt_covg,
                                                      gt = gt)

GenotypedVariant.objects.insert(gvs.values())





