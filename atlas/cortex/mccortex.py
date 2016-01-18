from __future__ import print_function

import os
import subprocess
import requests
import json
from Bio.Seq import Seq

class McCortexRunner(object):

  def __init__(self, sample, panels, seq, db_name, kmer = 31, force = False, panel_name = None):
    self.sample = sample
    self.panels = panels
    self.seq = seq
    self.db_name = db_name
    self.kmer = kmer      
    self.force = force
    self._panel_name = panel_name

  def run(self):
      if self.force or not os.path.exists(self.covg_tmp_file_path):
          self._check_panels() 
          self._run_cortex() 

  def _check_panels(self):
      ## If panel does not exists then build it
      for panel in self.panels:
          if not os.path.exists(panel.filepath):
              raise ValueError("Could not find a panel at %s." % panel.filepath)        

  def _run_cortex(self):
      ## If ctx binary does not exist then build it
      self._build_panel_binary_if_required()
      ## Now get coverage on panel
      self._run_coverage_if_required()

  def _build_panel_binary_if_required(self):
      if not os.path.exists(self.ctx_skeleton_filepath) or self.force:
          if os.path.exists(self.ctx_skeleton_filepath):
              os.remove(self.ctx_skeleton_filepath)
          ## panel  
          seq_list = self._create_sequence_list() 
          cmd = ["/home/phelimb/git/mccortex/bin/mccortex31", "build", "-q",
                   "-k", str(self.kmer)] + seq_list + [self.ctx_skeleton_filepath]
          # print (cmd)
          subprocess.check_output(cmd)
  def _create_sequence_list(self):
      seq_list = []
      for panel in self.panels:
         seq_list.extend(["-s", "%s" % panel.name, "-1", panel.filepath])
      return seq_list


  def _run_coverage_if_required(self):
      if not os.path.exists(self.ctx_tmp_filepath) or not os.path.exists(self.covg_tmp_file_path) or self.force:
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
             "-k", str(self.kmer), "-s", self.sample_name,
             "-o", self.covg_tmp_file_path]
      for seq in self.seq:
          cmd.extend(["-1", seq])
      for panel in self.panels:
      	  cmd.extend(["-c", panel.filepath])
      cmd.append(self.ctx_tmp_filepath)
      # print (cmd)
      return cmd 

  @property 
  def sample_name(self):
      return "-".join([self.sample, self.db_name, str(self.kmer)])

  @property
  def panel_name(self):
      if self._panel_name is None:
          self._panel_name = "-".join([p.name for p in self.panels])
      return self._panel_name

  @property 
  def sample_panel_name(self):
      return "_".join([self.sample_name, self.panel_name])

  @property
  def ctx_tmp_filepath(self):
      return "/tmp/%s.ctx" % self.sample_panel_name

  @property
  def covg_tmp_file_path(self):
      return "/tmp/%s.covgs" % self.sample_panel_name

  @property 
  def ctx_skeleton_filepath(self):
      return os.path.abspath("data/skeletons/%s_%i.ctx" % (self.panel_name, self.kmer))

class McCortexQuery(object):

  def __init__(self, port, seeds = []):
      self.base_url = "http://localhost:%i/" % port
      self.seeds = seeds

  def query(self, kmer):
      return McCortexQueryResult(kmer, requests.get(self.base_url + kmer).json(), seeds = self.seeds)

class McCortexQueryResult(object):


  def __init__(self, query_kmer, data, seeds = []):
      self.query_kmer = query_kmer
      self.data = data
      self.seeds = seeds

  def __str__(self):
      return str(self.data)

  def __repr__(self):
      return str(self.data)

  def forward(self):
      forward_kmers = []
      if self.complement:
          for l in self.left:
            forward_kmers.append(str(Seq(l + self.data["key"][:-1]).reverse_complement()))
      else:
          for r in self.right:
            forward_kmers.append(str(self.data["key"][1:] + r))
      if len(forward_kmers) > 1:
          forward_kmers = [f for f in forward_kmers if f in self.seeds]
      return forward_kmers

  @property
  def right(self):
      return self.data["right"]

  @property
  def left(self):
      return self.data["left"]     

  @property
  def kmer(self):
      return self.data["key"]

  @property
  def depth(self):
      try:
          return int(self.data["edges"])
      except ValueError:
          try:
              return int(self.data["edges"][-1])
          except ValueError:
              return 0
  @property
  def complement(self):
      return self.kmer != self.query_kmer

        

