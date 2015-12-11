from math import exp 
from math import factorial 
from math import log 

def percent_coverage_from_expected_coverage(coverage):
	## With low coverage we expect a lower percent of the sequence to be coverage. 
	return 1 - exp(-coverage)

def log_poisson_prob(lam, k):
	return -lam  + k*log(lam) - log_factorial(k); 

def log_factorial(n):
	assert n >= 0 
	out = 0
	for i in range(n):
		out += log(i + 1)
	return out

def log_lik_depth(depth, expected_depth):
	if expected_depth <= 0:
		raise ValueError("Expected depth must be greater than 0")
	if depth < 0:
		raise ValueError("Depth must not be negative")		
	return log_poisson_prob(lam = expected_depth, k = depth)

class Typer(object):

	def __init__(self):
		pass

	def type(self, l):
		raise NotImplemented("Implemented in sub class")

class GeneTyper(object):

	def __init__(self, depths, contamination_depths = []):
		self.depths = depths
		self.contamination_depths = contamination_depths

	def type(self, genes):
		typed_as_present = []
		for gene_name, gene_versions in genes.iteritems():
			best_gene_version = self.get_best_gene_version(gene_versions.values())
			if best_gene_version.percent_coverage > 30:
				typed_as_present.append(best_gene_version)
		return typed_as_present

	def get_best_gene_version(self, genes):
		genes.sort(key=lambda x: x.percent_coverage, reverse=True)
		current_best_gene = genes[0]
		for gene in genes[1:]:
			if gene.percent_coverage < current_best_gene.percent_coverage:
				return current_best_gene
			else:
				if gene.depth > current_best_gene.depth:
					current_best_gene = gene
		return current_best_gene