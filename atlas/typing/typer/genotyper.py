from __future__ import print_function
import os
import json
import csv
import glob
import logging
import subprocess

from atlas.schema import Variant
from atlas.typing import SequenceProbeCoverage
from atlas.typing import VariantProbeCoverage
from atlas.typing import ProbeCoverage
from atlas.typing import Panel

from atlas.typing.typer.presence import GeneCollectionTyper
from atlas.typing.typer.variant import VariantTyper

from atlas.panelgeneration import VariantPanel
from atlas.schema import VariantCallSet

from atlas.cortex import McCortexRunner

from atlas.utils import get_params

LOGGER = logging.getLogger("logger")


def max_pnz_threshold(vp):
    t = max(
        100 - 2 * math.floor(float(max([len(alt) for alt in vp.alts])) / 100), 30)
    return t


class CoverageParser(object):

    def __init__(
            self,
            args,
            panel_file_paths,
            panels=None,
            verbose=True,
            skeleton_dir='/tmp/'):
        self.args = args
        self.covgs = {"variant": {}, "presence": {}}
        self.variant_covgs = self.covgs["variant"]
        self.gene_presence_covgs = self.covgs["presence"]
        self.mc_cortex_runner = None
        self.verbose = verbose
        self.skeleton_dir = skeleton_dir
        self.panel_file_paths = panel_file_paths
        self.panels = []
        for panel_file_path in self.panel_file_paths:
            panel = Panel(panel_file_path)
            self.panels.append(panel)

    def run(self):
        self._run_cortex()
        self._parse_covgs()

    def _run_cortex(self):
        self.mc_cortex_runner = McCortexRunner(sample=self.args.sample,
                                               panels=self.panels,
                                               seq=self.args.seq,
                                               db_name=self.args.db_name,
                                               kmer=self.args.kmer,
                                               force=self.args.force,
                                               panel_name=self.panel_name,
                                               skeleton_dir=self.skeleton_dir)
        self.mc_cortex_runner.run()

    @property
    def panel_name(self):
        return "-".join([panel.name for panel in self.panels])

    def _parse_summary_covgs_row(self, row):
        try:
            return row[0], int(row[2]), int(row[3]), 100 * float(row[4])
        except ValueError:
            LOGGER.warning("Failed to parse %s" % str(row))
            return row[0], 0, 0, 0.0

    def _parse_covgs(self):
        with open(self.mc_cortex_runner.covg_tmp_file_path, 'r') as infile:
            self.reader = csv.reader(infile, delimiter="\t")
            for row in self.reader:
                allele, median_depth, min_depth, percent_coverage = self._parse_summary_covgs_row(
                    row)
                allele_name = allele.split('?')[0]
                if self._is_variant_panel(allele_name):
                    self._parse_variant_panel(row)
                else:
                    self._parse_seq_panel(row)

    def _is_variant_panel(self, allele_name):
        try:
            alt_or_ref, _id = allele_name.split('-')
            return bool(alt_or_ref)
        except ValueError:
            return False

    def _parse_seq_panel(self, row):
        allele, median_depth, min_depth, percent_coverage = self._parse_summary_covgs_row(
            row)
        allele_name = allele.split('?')[0]
        params = get_params(allele)
        panel_type = params.get("panel_type", "presence")
        name = params.get('name')
        if panel_type in ["variant", "presence"]:
            gp = SequenceCoverage.create_object(
                name=name,
                version=params.get(
                    'version',
                    'NA'),
                percent_coverage=percent_coverage,
                median_depth=median_depth,
                min_depth=min_depth,
                length=params.get("length"))
            try:
                self.covgs[panel_type][gp.name][gp.version] = gp
            except KeyError:
                self.covgs[panel_type][gp.name] = {}
                self.covgs[panel_type][gp.name][gp.version] = gp

        else:
            # Species panels are treated differently
            l = int(params.get("length", -1))
            try:
                self.covgs[panel_type][name]["total_bases"] += l
                if percent_coverage > 75 and median_depth > 0:
                    self.covgs[panel_type][name][
                        "percent_coverage"].append(percent_coverage)
                    self.covgs[panel_type][name]["length"].append(l)
                    self.covgs[panel_type][name]["median"].append(median_depth)
            except KeyError:
                if panel_type not in self.covgs:
                    self.covgs[panel_type] = {}
                self.covgs[panel_type][name] = {}
                self.covgs[panel_type][name]["total_bases"] = l
                if percent_coverage > 75 and median_depth > 0:
                    self.covgs[panel_type][name][
                        "percent_coverage"] = [percent_coverage]
                    self.covgs[panel_type][name]["length"] = [l]
                    self.covgs[panel_type][name]["median"] = [median_depth]
                else:
                    self.covgs[panel_type][name]["percent_coverage"] = []
                    self.covgs[panel_type][name]["length"] = []
                    self.covgs[panel_type][name]["median"] = []

    def _parse_variant_panel(self, row):
        allele, reference_median_depth, min_depth, reference_percent_coverage = self._parse_summary_covgs_row(
            row)
        var_name = allele.split('?')[0].split('-')[1]
        params = get_params(allele)
        num_alts = int(params.get("num_alts", 0))
        reference_coverage = ProbeCoverage(
            percent_coverage=reference_percent_coverage,
            median_depth=reference_median_depth,
            min_depth=min_depth)
        alternate_coverages = []
        for i in range(num_alts):
            row = self.reader.next()
            alt_allele, alternate_median_depth, min_depth, alternate_percent_coverage = self._parse_summary_covgs_row(
                row)
            alternate_coverages.append(
                ProbeCoverage(
                    min_depth=min_depth,
                    percent_coverage=alternate_percent_coverage,
                    median_depth=alternate_median_depth))
        variant_probe_coverage = VariantProbeCoverage(
            reference_coverage=reference_coverage,
            alternate_coverages=alternate_coverages,
            var_name=var_name,
            params=params)
        try:
            self.variant_covgs[allele].append(variant_probe_coverage)
        except KeyError:
            self.variant_covgs[allele] = [variant_probe_coverage]


class Genotyper(object):

    """Takes output of mccortex coverages and types"""

    def __init__(
            self,
            args,
            expected_depths,
            variant_covgs,
            gene_presence_covgs,
            contamination_depths=[],
            verbose=False,
            base_json={}):
        self.args = args
        self.variant_covgs = variant_covgs
        self.gene_presence_covgs = gene_presence_covgs
        self.out_json = base_json
        self.verbose = verbose
        self.expected_depths = expected_depths
        self.contamination_depths = contamination_depths

    def run(self):
        self._type()
        if not self.args.quiet:
            print(json.dumps(self.out_json,
                             indent=4, separators=(',', ': ')))

    def _type(self):
        self._type_genes()
        self._type_variants()

    def _type_genes(self):
        gt = GeneCollectionTyper(
            expected_depths=self.expected_depths,
            contamination_depths=self.contamination_depths)
        gene_presence_covgs_out = {}
        for gene_name, gene_collection in self.gene_presence_covgs.items():
            self.gene_presence_covgs[gene_name] = gt.genotype(gene_collection)
            gene_presence_covgs_out[
                gene_name] = self.gene_presence_covgs[gene_name].to_dict()
        self.out_json[self.args.sample][
            "typed_presence"] = gene_presence_covgs_out

    def _type_variants(self):
        self.out_json[self.args.sample]["typed_variants"] = {}
        out_json = self.out_json[self.args.sample]["typed_variants"]
        gt = VariantTyper(
            expected_depths=self.expected_depths,
            contamination_depths=self.contamination_depths)

        for variant, probe_coverages in self.variant_covgs.items():
            call = gt.type(probe_coverages)
            if sum(call.genotype) > 0:
                out_json[variant] = call.to_mongo().to_dict()
        # typed_variants = gt.type(self.variant_covgs)
        # self.out_json[self.args.sample]["typed_variants"] = {}

        # for name, tvs in typed_variants.items():
        #     for tv in tvs:
        #         try:
        #             out_json[name].append(tv.to_dict())
        #         except KeyError:
        #             out_json[name] = [tv.to_dict()]
