from atlas.typing.typer.base import Typer

class VariantTyper(Typer):

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