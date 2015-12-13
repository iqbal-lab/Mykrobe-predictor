#! /usr/bin/env python
import sys
import os
import json
from pprint import pprint
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
from atlas.vcf2db import TypedVariant
from atlas.vcf2db import VariantPanel
from atlas.utils import check_args

from atlas.genotyping import GeneTyper
from atlas.genotyping import Presence as GenePresence

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

     

  @property 
  def gene_version(self):
      return "-".join([self.name, self.version])

  def to_dict(self):
      d  = {  "name" : self.name,
              "alt_name" : ",".join(self.alt_names),
              "covg" : {"percent_coverage" : self.percent_coverage, 
                        "median_depth" : self.depth
                        }
            }
      return d      

def build_binaries(args):
    ## If panel does not exists then build it
    panel_filepath = os.path.abspath("data/panels/%s.fasta" % args.panel) 
    if not os.path.exists(panel_filepath):
        raise ValueError("Could not find a panel at %s. Run 'atlas dump'. " % panel_filepath)
    ## If ctx binary does not exist then build it
    ctx_skeleton_filepath = os.path.abspath("data/skeletons/%s.ctx" % args.panel) 
    if not os.path.exists(ctx_skeleton_filepath):
        subprocess.check_output(["/home/phelimb/git/mccortex/bin/mccortex31", "build", "-q",
                                 "-k", str(args.kmer), "-s", "%s" % args.panel,
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
    else:
        # print "Warning: Using pre-built binaries. Run with --force if panel has been updated."
        pass  

def run(parser, args):
    args = parser.parse_args()
    args = check_args(args)
    connect('atlas-%s-%i' % (args.db_name ,args.kmer))
    if not args.panel:
        args.panel = "panel-%s-%i" % (args.db_name, args.kmer)
    try:
        call_set = CallSet.objects.get(name = args.sample + "_%s" % args.name)
    except DoesNotExist:
        call_set = CallSet.create(name = args.sample  + "_%s" % args.name, sample_id = args.sample)
    ## Clear any genotyped calls so far
    TypedVariant.objects(call_set = call_set).delete()
    gvs = {}
    gene_presence = {}

    with open("/tmp/%s.covg" % sample, 'r') as infile:
        reader = csv.reader(infile, delimiter = "\t")
        for row in reader:
            allele = row[0]
            median_coverage = int(row[2])
            pnz = 100 * float(row[3])
            allele_name = allele.split('?')[0]
            try:
                alt_or_ref, _id = allele_name.split('-')
            except ValueError: 
                alt_or_ref = None
                _id = allele_name
            # vp = VariantPanel.objects.get(id = _id)
            MAX_PNZ_THRESHOLD = 1.0#max_pnz_threshold(vp)
            if not alt_or_ref:
                if pnz > 0:
                    params = get_params(allele)
                    gp = GenePresence(name = params.get('name'),
                                 version = params.get('version', 'N/A'),
                                 percent_coverage = pnz,
                                 depth = median_coverage
                                 )
                    try:
                        gene_presence[gp.name][gp.version] = gp
                    except KeyError:
                        gene_presence[gp.name] = {}
                        gene_presence[gp.name][gp.version] = gp
            else:               
                if alt_or_ref == "ref":
                    params = get_params(allele)
                    ref_pnz = pnz
                    ref_covg = median_coverage

                    num_alts = int(params.get("num_alts"))
                    alt_pnz = 0
                    alt_covg = 0            
                    for _ in range(num_alts):
                      row = reader.next()
                      allele = row[0]
                      median_coverage = int(row[2])
                      pnz = 100 * float(row[3])      
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
                else:
                    gt = "-/-"
                    # print(alt_pnz, ref_pnz, alt_covg, ref_covg)
                if gt not in  ["0/0", "-/-"] or args.all:
                   # print _id, params.get("gene"), params.get("mut"), ref_covg, alt_covg, ref_pnz, alt_pnz, gt
                   gvs[_id] = TypedVariant.create_object(name = _id,
                                                              call_set = call_set,
                                                              ref_pnz = ref_pnz, 
                                                              alt_pnz = alt_pnz,
                                                              ref_coverage = ref_covg, 
                                                              alt_coverage = alt_covg,
                                                              gt = gt,
                                                              alt_name = "_".join([params.get("gene"), params.get("mut")]))
    # print(json.dumps({args.sample : [gv.to_dict() for gv in gvs.values()]},
    #                   indent=4, separators=(',', ': ')))
    gt = GeneTyper(depths = [100])
    gt.type(gene_presence)
    print(json.dumps({args.sample : [gv.to_dict() for gv in gt.type(gene_presence)]},
                      indent=4, separators=(',', ': ')))    

    # print(json.dumps({args.sample : {"gene_presence" : gene_presngene_presence}},
    #                   indent=4, separators=(',', ': ')))

    # TypedVariant.objects.insert(gvs.values())