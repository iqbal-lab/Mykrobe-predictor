from __future__ import print_function
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
		# self._load_taxon_thresholds()
		self._aggregate_all()

	def _load_taxon_thresholds(self):
		with open("data/predict/taxon_coverage_threshold.json", "r") as infile:
			self.threshold = json.load(infile)

	def _aggregate_all(self):
		self.phylo_group_summary = self._aggregate(self.phylo_group_covgs)
		self.species_summary = self._aggregate(self.species_covgs)
		self.lineage_summary = self._aggregate(self.lineage_covgs)
		self.out_json["phylogenetics"] = {}
		self.out_json["phylogenetics"]["phylo"] = self.phylo_group_covgs
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

