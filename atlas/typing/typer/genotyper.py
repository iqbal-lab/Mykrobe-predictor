import datetime
import os
import json
from pprint import pprint
import csv
import glob
from mongoengine import connect
from mongoengine import DoesNotExist
import subprocess

from atlas.typing import TypedVariant
from atlas.typing import TypedPresence
from presence import PresenceTyper
from variant import VariantTyper
from atlas.vcf2db import VariantPanel
from atlas.vcf2db import CallSet

def get_params(url):
    params = {}
    p_str = url.split("?")[1]
    p_str = p_str.split('&')
    for p in p_str:
        k,v = p.split("=")
        params[k] = v
    return params

def max_pnz_threshold(vp):
    t =  max(100 - 2 * math.floor(float(max([len(alt) for alt in vp.alts])) / 100), 30)
    return t

class Genotyper(object):

  def __init__(self, args):
    self.args = args
    self.gvs = {}
    self.gene_presence = {}
    if not args.panel:
        args.panel = "panel-%s-%i" % (args.db_name, args.kmer)

  def run(self):
      self._check_args() 
      self._run_cortex() 
      self._parse_covgs()
      self._type()      
      self._connect_to_db()
      self._set_up_db()
      # self._insert_to_db()


  def _type(self):
      self._type_genes()
      self._type_variants()

  def _type_genes(self):
      gt = PresenceTyper(depths = [100])
      gt.type(self.gene_presence)
      print(json.dumps({self.args.sample : [gv.to_dict() for gv in gt.type(self.gene_presence)]},
                        indent=4, separators=(',', ': '))) 

  def _type_variants(self):
      pass 

  def _parse_covgs(self):
    with open(self.covg_tmp_file_path, 'r') as infile:
        reader = csv.reader(infile, delimiter = "\t")
        for row in reader:
            allele, median_depth, percent_coverage = row[0], int(row[2]), 100*float(row[3])
            allele_name = allele.split('?')[0]
            try:
                alt_or_ref, _id = allele_name.split('-')
            except ValueError: 
                alt_or_ref = None
                _id = allele_name
            # vp = VariantPanel.objects.get(id = _id)
            MAX_PNZ_THRESHOLD = 1.0#max_pnz_threshold(vp)
            if not alt_or_ref:
                if percent_coverage > 0:
                    params = get_params(allele)
                    gp = TypedPresence(name = params.get('name'),
                                 version = params.get('version', 'N/A'),
                                 percent_coverage = percent_coverage,
                                 median_depth = median_depth
                                 )
                    try:
                        self.gene_presence[gp.name][gp.version] = gp
                    except KeyError:
                        self.gene_presence[gp.name] = {}
                        self.gene_presence[gp.name][gp.version] = gp
            else:               
                if alt_or_ref == "ref":
                    params = get_params(allele)
                    reference_percent_coverage = pnz
                    reference_median_depth = median_depth

                    num_alts = int(params.get("num_alts"))
                    alternate_percent_coverage = 0
                    alternate_median_depth = 0            
                    for _ in range(num_alts):
                      row = reader.next()
                      allele = row[0]
                      median_depth = int(row[2])
                      pnz = 100 * float(row[3])      
                      if pnz >= alternate_percent_coverage:
                          alternate_percent_coverage = pnz
                          if median_depth > alternate_median_depth:
                              alternate_median_depth = median_depth          
                if alternate_percent_coverage >= MAX_PNZ_THRESHOLD and reference_percent_coverage < MAX_PNZ_THRESHOLD and alternate_median_depth > 0:
                    gt = "1/1"
                elif alternate_percent_coverage >= MAX_PNZ_THRESHOLD and reference_percent_coverage >= MAX_PNZ_THRESHOLD and alternate_median_depth > 0:
                    gt = "0/1"
                elif alternate_median_depth < MAX_PNZ_THRESHOLD:
                    gt = "0/0"
                else:
                    gt = "-/-"
                    # print(alternate_percent_coverage, reference_percent_coverage, alternate_median_depth, reference_median_depth)
                if gt not in  ["0/0", "-/-"] or self.args.all:
                   # print _id, params.get("gene"), params.get("mut"), reference_median_depth, alternate_median_depth, reference_percent_coverage, alternate_percent_coverage, gt
                   self.gvs[_id] = TypedVariant.create_object(name = _id,
                                                              call_set = self.call_set,
                                                              reference_percent_coverage = reference_percent_coverage, 
                                                              alternate_percent_coverage = alternate_percent_coverage,
                                                              ref_coverage = reference_median_depth, 
                                                              alt_coverage = alternate_median_depth,
                                                              gt = gt,
                                                              alt_name = "_".join([params.get("gene"), params.get("mut")]))  
  def _connect_to_db(self):
    connect('atlas-%s-%i' % (self.args.db_name ,self.args.kmer))

  def _set_up_db(self):
    try:
        self.call_set = CallSet.objects.get(name = self.args.sample + "_%s" % self.args.name)
    except DoesNotExist:
        self.call_set = CallSet.create(name = self.args.sample  + "_%s" % self.args.name, sample_id = self.args.sample)
    ## Clear any genotyped calls so far
    TypedVariant.objects(call_set = self.call_set).delete()

  def _build_panel_binary_if_required(self):
      if not os.path.exists(self.ctx_skeleton_filepath) or self.args.force:
          if os.path.exists(self.ctx_skeleton_filepath):
              os.remove(self.ctx_skeleton_filepath)        
          subprocess.check_output(["/home/phelimb/git/mccortex/bin/mccortex31", "build", "-q",
                                   "-k", str(self.args.kmer), "-s", "%s" % self.args.panel,
                                   "-1", self.panel_filepath, self.ctx_skeleton_filepath])    

  def _check_args(self):
      ## If panel does not exists then build it
      if not os.path.exists(self.panel_filepath):
          raise ValueError("Could not find a panel at %s. Run 'atlas dump'. " % self.panel_filepath)    

  def _run_cortex(self):
      ## If ctx binary does not exist then build it
      self._build_panel_binary_if_required()
      ## Now get coverage on panel
      self._run_coverage_if_required()

  def _run_coverage_if_required(self):
      if not os.path.exists(self.ctx_tmp_filepath) or not os.path.exists(self.covg_tmp_file_path) or self.args.force:
          if os.path.exists(self.ctx_tmp_filepath):
              os.remove(self.ctx_tmp_filepath)
          if os.path.exists(self.covg_tmp_file_path):
              os.remove(self.covg_tmp_file_path)      
          subprocess.check_output(self.coverages_cmd)
      else:
          # print "Warning: Using pre-built binaries. Run with --force if panel has been updated."
          pass  
  
  @property 
  def coverages_cmd(self):
      cmd = ["/home/phelimb/git/mccortex/bin/mccortex31", "geno", "-q",
             "-I", self.ctx_skeleton_filepath,
             "-k", str(self.args.kmer), "-s", self.sample,
             "-o", self.covg_tmp_file_path]
      for seq in self.args.seq:
        cmd.extend(["-1", seq])
      cmd.extend(["-c", self.panel_filepath, self.ctx_tmp_filepath])    
      return cmd

  @property 
  def sample(self):
      return "-".join([self.args.sample, self.args.db_name, str(self.args.kmer)])

  @property 
  def ctx_skeleton_filepath(self):
    return os.path.abspath("data/skeletons/%s_%i.ctx" % (self.args.panel, self.args.kmer)) 

  @property
  def ctx_tmp_filepath(self):
    return "/tmp/%s_%s.ctx" % (self.sample, self.args.panel)

  @property
  def covg_tmp_file_path(self):
      return "/tmp/%s_%s.covgs" % (self.sample, self.args.panel)

  @property
  def panel_filepath(self):
      return os.path.abspath("data/panels/%s.fasta" % self.args.panel) 