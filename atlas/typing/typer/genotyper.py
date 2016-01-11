from __future__ import print_function
import os
import json
from pprint import pprint
import csv
import glob
from mongoengine import connect
from mongoengine import DoesNotExist
import subprocess

from atlas.typing import TypedVariant
from atlas.typing import SequenceCoverage
from atlas.typing import Panel

from atlas.typing.typer.presence import GeneCollectionTyper
from atlas.typing.typer.variant import VariantTyper

from atlas.vcf2db import VariantPanel
from atlas.vcf2db import CallSet

from atlas.cortex import McCortexRunner


def get_params(url):
    params = {}
    try:
        p_str = url.split("?")[1]
    except IndexError:
        return params
    p_str = p_str.split('&')
    for p in p_str:
        k,v = p.split("=")
        params[k] = v
    return params

def max_pnz_threshold(vp):
    t =  max(100 - 2 * math.floor(float(max([len(alt) for alt in vp.alts])) / 100), 30)
    return t

class CoverageParser(object):

  def __init__(self, args, panels = None, verbose = True):
    self.args = args
    self.covgs = {"variant" : {}, "presence" : {}}
    self.variant_covgs = self.covgs["variant"]
    self.gene_presence_covgs = self.covgs["presence"]
    self.out_json = {self.args.sample : {}}   
    self.mc_cortex_runner = None
    self.verbose = verbose
    if panels:
        self.panel_names = panels
    else:
        self.panel_names = ["panel-%s-%i" % (args.db_name, args.kmer)] 

  def run(self):
      self._connect_to_db()      
      self._set_up_db()
      self._run_cortex()
      self._parse_covgs()       

  def _connect_to_db(self):
    connect('atlas-%s-%i' % (self.args.db_name ,self.args.kmer))

  def _set_up_db(self):
    try:
        self.call_set = CallSet.objects.get(name = self.args.sample + "_%s" % self.args.name)
    except DoesNotExist:
        self.call_set = CallSet.create(name = self.args.sample  + "_%s" % self.args.name, sample_id = self.args.sample)
    ## Clear any genotyped calls so far
    TypedVariant.objects(call_set = self.call_set).delete()  

  def _run_cortex(self):
      self.mc_cortex_runner = McCortexRunner(sample = self.args.sample,
                                             panels = self.panels,
                                             seq = self.args.seq,
                                             db_name = self.args.db_name,
                                             kmer = self.args.kmer,
                                             force = self.args.force)
      self.mc_cortex_runner.run()

  @property
  def panels(self):
      panels = []
      for panel in self.panel_names:
          panels.append(Panel(panel))
      return panels

  def _parse_summary_covgs_row(self, row):
      return row[0], int(row[2]), 100*float(row[3])

  def _parse_covgs(self):
      with open(self.mc_cortex_runner.covg_tmp_file_path, 'r') as infile:
          self.reader = csv.reader(infile, delimiter = "\t")
          for row in self.reader:
              allele, median_depth, percent_coverage = self._parse_summary_covgs_row(row)
              allele_name = allele.split('?')[0]
              if self._is_variant_panel(allele_name):
                  self._parse_variant_panel(row)
              else:
                  self._parse_seq_panel(row)

  def _is_variant_panel(self, allele_name):
      try:
        alt_or_ref, _id = allele_name.split('-')
        return bool(alt_or_ref)      
      except ValueError, e:
        return False


  def _parse_seq_panel(self, row):
      allele, median_depth, percent_coverage = self._parse_summary_covgs_row(row)
      allele_name = allele.split('?')[0]    
      params = get_params(allele)
      panel_type = params.get("panel_type", "presence")
      gp = SequenceCoverage.create_object(name = params.get('name'),
                   version = params.get('version', 'N/A'),
                   percent_coverage = percent_coverage,
                   median_depth = median_depth,
                   length = params.get("length")
                   )
      try:
          self.covgs[panel_type][gp.name][gp.version] = gp
      except KeyError:
          try:
              self.covgs[panel_type][gp.name] = {}
          except KeyError:
              self.covgs[panel_type] = {}
              self.covgs[panel_type][gp.name] = {}
          finally:
              self.covgs[panel_type][gp.name][gp.version] = gp

  def _parse_variant_panel(self, row):
      allele, reference_median_depth, reference_percent_coverage = self._parse_summary_covgs_row(row)
      allele_name = allele.split('?')[0].split('-')[1]
      params = get_params(allele)   
      num_alts = int(params.get("num_alts"))
      for i in range(num_alts):
          row = self.reader.next()
          allele, alternate_median_depth, alternate_percent_coverage = self._parse_summary_covgs_row(row)
          try:
              alt_name = "_".join([params["gene"], params["mut"]])
          except KeyError:
              alt_name = ""
          tv = TypedVariant.create_object(
                                      name = allele_name,
                                      call_set = self.call_set,
                                      reference_percent_coverage = reference_percent_coverage, 
                                      alternate_percent_coverage = alternate_percent_coverage,
                                      reference_median_depth = reference_median_depth, 
                                      alternate_median_depth = alternate_median_depth,
                                      alt_name = alt_name,
                                      alt_index = i)
          try:
              self.variant_covgs[allele_name].append(tv)
          except KeyError:
              self.variant_covgs[allele_name] = [tv]



class Genotyper(object):

  """Takes output of mccortex coverages and types"""

  def __init__(self, args, depths, variant_covgs, gene_presence_covgs, verbose = False):
    self.args = args
    self.variant_covgs = variant_covgs
    self.gene_presence_covgs = gene_presence_covgs
    self.out_json = {self.args.sample : {}}   
    self.verbose = verbose
    self.depths = depths

  def run(self):
      self._type()    
      if not self.args.quiet:
          print(json.dumps(self.out_json,
                        indent=4, separators=(',', ': ')))           
      # self._insert_to_db()

  def _type(self):
      self._type_genes()
      self._type_variants()

  def _type_genes(self):
      gt = GeneCollectionTyper(depths = self.depths)
      gene_presence_covgs_out = {}
      for gene_name, gene_collection in self.gene_presence_covgs.iteritems():
          self.gene_presence_covgs[gene_name] = gt.genotype(gene_collection)
          if self.verbose or self.gene_presence_covgs[gene_name].gt not in ["0/0", "-/-"]:
              gene_presence_covgs_out[gene_name] = self.gene_presence_covgs[gene_name].to_dict()
      self.out_json[self.args.sample]["typed_presence"]  = gene_presence_covgs_out

  def _type_variants(self):
      gt = VariantTyper(depths = self.depths) 
      typed_variants = gt.type(self.variant_covgs)
      self.out_json[self.args.sample]["typed_variants"] = {}
      out_json = self.out_json[self.args.sample]["typed_variants"] 
      for name, tvs in typed_variants.iteritems():
          for tv in tvs:
              if self.verbose or tv.gt not in ["0/0", "-/-"]:
                  try:
                      out_json[name].append(tv.to_dict())
                  except KeyError:
                      out_json[name] = [tv.to_dict()]

                                             

