from atlas.typing.typer.base import Typer
from atlas.stats import log_lik_depth
from atlas.stats import percent_coverage_from_expected_coverage
from math import log
MIN_LLK = -99999999


class PresenceTyper(Typer):

    "Initiated with expected depths and contamination depths"

    def __init__(self, depths, contamination_depths = []):
        super(PresenceTyper, self).__init__(depths, contamination_depths)

    def genotype(self, sequence_coverage):
        "Takes a single SequenceCoverage object (or child) and returns genotype"
        if not self.has_contamination():
            self._type_with_minor_model(sequence_coverage)
        else:
            self._type_without_minor_model(sequence_coverage)
        return sequence_coverage

    def _type_with_minor_model(self, sequence_coverage):
        for expected_depth in self.depths:
            hom_alt_likelihood = log_lik_depth(sequence_coverage.median_depth,
                                               expected_depth * 0.75)
            het_likelihood     = log_lik_depth(sequence_coverage.median_depth,
                                               expected_depth * self.minimum_detectable_frequency )
            hom_ref_likelihood = log_lik_depth(sequence_coverage.median_depth, expected_depth * 0.001)
            ## Posterior
            hom_ref_likelihood = self._log_post_hom_ref(hom_ref_likelihood)
            hom_alt_likelihood = self._log_post_het_or_alt(hom_alt_likelihood,
            											   expected_depth * 0.75,
            											   sequence_coverage)
            het_likelihood = self._log_post_het_or_alt(het_likelihood,
            										  expected_depth * self.minimum_detectable_frequency,
            										  sequence_coverage)

        gt = self.likelihoods_to_genotype([hom_ref_likelihood,
                                           het_likelihood,
                                           hom_alt_likelihood])
        sequence_coverage.set_genotype(gt)

    def _type_without_minor_model(self, sequence_coverage):
        if sequence_coverage.percent_coverage > sequence_coverage.percent_coverage_threshold:
            sequence_coverage.set_genotype("1/1")
        else:
            sequence_coverage.set_genotype("0/0")

    @property
    def minimum_detectable_frequency(self):
        if self.error_rate < 0.1:
            return 0.05
        else:
            return 0.25

    def _log_post_hom_ref(self, llk):
        return log(1) + llk

    def _log_post_het_or_alt(self, llk, expected_depth, sequence_coverage):
        expected_percentage_coverage = percent_coverage_from_expected_coverage(expected_depth)
        minimum_percentage_coverage_required =  expected_percentage_coverage * sequence_coverage.percent_coverage_threshold
        if sequence_coverage.percent_coverage > minimum_percentage_coverage_required:
            return self._log_post_hom_ref(llk)
        else:
            return MIN_LLK


class GeneCollectionTyper(Typer):

    """Types a collection of genes returning only the most likely version 
        in the collection"""

    def __init__(self, depths, contamination_depths = []):
        super(GeneCollectionTyper, self).__init__(depths, contamination_depths)
        self.presence_typer = PresenceTyper(depths, contamination_depths)

    def genotype(self, sequence_coverage_collection):
        """Types a collection of genes returning the most likely gene version 
            in the collection with it's genotype"""
        best_version = self.get_best_version(sequence_coverage_collection.values())
        return self.presence_typer.genotype(best_version)

    def get_best_version(self, sequence_coverages):
        sequences.sort(key=lambda x: x.percent_coverage, reverse=True)
        current_best_gene = sequences[0]
        for gene in sequences[1:]:
            if gene.percent_coverage < current_best_gene.percent_coverage:
                return current_best_gene
            else:
                if gene.median_depth > current_best_gene.median_depth:
                    current_best_gene = gene
        return current_best_gene

