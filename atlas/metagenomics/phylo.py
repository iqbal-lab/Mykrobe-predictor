from __future__ import print_function
import json
DEFAULT_THRESHOLD = 30

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

class SpeciesPredictor(object):

    def __init__(self, phylo_group_covgs, species_covgs, lineage_covgs, base_json):
        self.phylo_group_covgs = phylo_group_covgs
        self.species_covgs = species_covgs
        self.lineage_covgs = lineage_covgs
        self.out_json = base_json
        self.threshold = {}

    def run(self):
        self._load_taxon_thresholds()
        self._aggregate_all()

    def _load_taxon_thresholds(self):
        with open("data/predict/taxon_coverage_threshold.json", "r") as infile:
            self.threshold = json.load(infile)

    def _aggregate_all(self):
        self._aggregate(self.phylo_group_covgs)
        self._aggregate(self.species_covgs)
        self._aggregate(self.lineage_covgs)
        self.out_json["phylogenetics"] = {}
        self.out_json["phylogenetics"]["phylo_group"] = self.phylo_group_covgs
        self.out_json["phylogenetics"]["species"] = self.species_covgs
        self.out_json["phylogenetics"]["lineage"] = self.lineage_covgs

    def _aggregate(self, covgs):
        del_nodes = []
        for node, covg_collection  in covgs.iteritems():
            bases_covered = covg_collection["bases_covered"]
            total_bases = covg_collection["total_bases"]
            _median = covg_collection["median"]            
            aggregate_percent_covg = bases_covered/total_bases
            if aggregate_percent_covg >= self.threshold.get(node, DEFAULT_THRESHOLD):
                covgs[node] = {"percent_coverage" : bases_covered/total_bases, "median_depth" : median(_median)}
            else:
                del_nodes.append(node)
        for node in del_nodes:
            del covgs[node]

class AMRSpeciesPredictor(SpeciesPredictor):

    def __init__(self, phylo_group_covgs, species_covgs, lineage_covgs,
                 base_json):
        super(AMRSpeciesPredictor, self).__init__(phylo_group_covgs, species_covgs, lineage_covgs,
                 base_json)

    def is_saureus_present(self):
        return "Staphaureus" in self.out_json["phylogenetics"]["phylo_group"]

    def is_mtbc_present(self):
        return "MTBC" in self.out_json["phylogenetics"]["phylo_group"]

    def is_gram_neg_present(self):
        return self.is_klebsiella_pneumoniae_present() or self.is_escherichia_coli_present()
                
    def is_klebsiella_pneumoniae_present(self):
        return "Klebsiella_pneumoniae" in self.out_json["phylogenetics"]["species"]

    def is_escherichia_coli_present(self):
        return "Escherichia_coli" in self.out_json["phylogenetics"]["species"]

    def contamination_depths(self):
        contamination_depths = []
        ignore = []
        if self.is_saureus_present():
            ignore.append("Saureus")
        elif self.is_mtbc_present():
            ignore.append("Mtuberculosis")
        elif self.is_escherichia_coli_present():
            ignore.append("Escherichia_coli")
        elif self.is_klebsiella_pneumoniae_present():
            ignore.append("Klebsiella_pneumoniae")
        for node, covg_collection in self.species_covgs.iteritems():
            if node not in ignore:
                contamination_depths.append(covg_collection["median_depth"])
        return contamination_depths






