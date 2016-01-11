from atlas.typing.typer.base import Typer
from atlas.stats import log_lik_R_S_coverage
from atlas.typing.typer.base import MIN_LLK


class VariantTyper(Typer):

	def __init__(self, depths, contamination_depths = [], error_rate = 0.05):
		super(VariantTyper, self).__init__(depths, contamination_depths, error_rate)
		self.method = "MAP"

	def type(self, variants):
		if not self.has_contamination():
			self._type_with_minor_model(variants)
		else:
			self._type_without_minor_model(variants)
		return variants

	def _type_with_minor_model(self, variants):
		for variant_name, variants in variants.iteritems():
			for variant in variants:
				self._type_individual(variant, include_minor = True)

		return variants

	def _type_individual(self, variant, include_minor = True):
		hom_ref_likelihood = self._hom_ref_lik(variant)
		hom_alt_likelihood = self._hom_alt_lik(variant)
		if include_minor:
			het_likelihood = self._het_lik(variant)
		else:
			het_likelihood = MIN_LLK
		gt = self.likelihoods_to_genotype([hom_ref_likelihood, het_likelihood, hom_alt_likelihood])
		variant.set_genotype(gt)
		variant.set_copy_number(float(variant.alternate_median_depth) / self.depths[0])		

	def _hom_ref_lik(self, variant):
		if variant.reference_percent_coverage < 100:
			return MIN_LLK
		else:			
			hom_ref_likes = []
			## Either alt+cov or alt_covg + contam_covg
			for expected_depth in self.depths:
				hom_ref_likes.append(log_lik_R_S_coverage(variant.reference_median_depth, 
														  variant.alternate_median_depth,
														  expected_depth,
														  error_rate = self.error_rate))			
				for contamination in self.contamination_depths:
					hom_ref_likes.append(log_lik_R_S_coverage(variant.reference_median_depth, 
															  variant.alternate_median_depth,
															  expected_depth + contamination,
															  error_rate = self.error_rate))				
			return max(hom_ref_likes)

	def _hom_alt_lik(self, variant):
		if variant.alternate_percent_coverage < 100:
			return MIN_LLK
		else:	
			hom_alt_liks = []
			## Either alt+cov or alt_covg + contam_covg
			for expected_depth in self.depths:
				hom_alt_liks.append(log_lik_R_S_coverage(variant.alternate_median_depth,
														  variant.reference_median_depth, 
														  expected_depth,
														  error_rate = self.error_rate))			
				for contamination in self.contamination_depths:
					hom_alt_liks.append(log_lik_R_S_coverage(variant.alternate_median_depth,
															  variant.reference_median_depth,   
															  expected_depth + contamination,
															  error_rate = self.error_rate))				
			return max(hom_alt_liks)

	def _het_lik(self, variant):
		if variant.alternate_percent_coverage < 100 or variant.reference_percent_coverage < 100:
			return MIN_LLK
		else:
			return log_lik_R_S_coverage(variant.alternate_median_depth,
									    variant.reference_median_depth,   
									    self.depths[0],
									    error_rate = 0.1)

	def _type_without_minor_model(self, variants):
		raise NotImplementedError("")


