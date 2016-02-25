from atlas.typing.typer.base import Typer
from atlas.typing.typer.base import MIN_LLK
from atlas.schema import SequenceCall
from atlas.stats import log_lik_depth
from atlas.stats import percent_coverage_from_expected_coverage
from math import log


class PresenceTyper(Typer):

    "Initiated with expected depths and contamination depths"

    def __init__(self, expected_depths, contamination_depths=[]):
        super(
            PresenceTyper,
            self).__init__(
            expected_depths,
            contamination_depths)
        if len(expected_depths) > 1:
            raise NotImplementedError("Mixed samples not supported")

    def type(self, sequence_probe_coverage):
        "Takes a single SequenceCoverage object (or child) and returns genotype"
        call = self._type(sequence_probe_coverage)
        return call

    def _type(self, sequence_probe_coverage):
        hom_alt_likelihoods = []
        het_likelihoods = []
        hom_ref_likelihoods = []
        for expected_depth in self.expected_depths:
            hom_alt_likelihoods.append(
                self._hom_alt_likeihood(
                    median_depth=sequence_probe_coverage.median_depth,
                    expected_depth=expected_depth))
            if not self.has_contamination():
                het_likelihoods.append(
                    self._het_likelihood(
                        median_depth=sequence_probe_coverage.median_depth,
                        expected_depth=expected_depth))
            else:
                het_likelihoods.append(MIN_LLK)

            hom_ref_likelihoods.append(
                self._hom_ref_likelihood(
                    median_depth=sequence_probe_coverage.median_depth,
                    expected_depth=expected_depth))

            for contamination_depth in self.contamination_depths:
                hom_alt_likelihoods.append(
                    self._hom_alt_likeihood(
                        median_depth=sequence_probe_coverage.median_depth,
                        expected_depth=expected_depth +
                        contamination_depth))
                # NOTE : _HOM_ALT_LIKEIHOOD is not a typo
                hom_ref_likelihoods.append(
                    self._hom_alt_likeihood(
                        median_depth=sequence_probe_coverage.median_depth,
                        expected_depth=contamination_depth))
            # Posterior
        hom_ref_likelihood = self._log_post_hom_ref(max(hom_ref_likelihoods))
        hom_alt_likelihood = self._log_post_het_or_alt(
            max(hom_alt_likelihoods),
            expected_depth * 0.75,
            sequence_probe_coverage)
        het_likelihood = self._log_post_het_or_alt(
            max(het_likelihoods),
            expected_depth *
            self.minimum_detectable_frequency,
            sequence_probe_coverage)
        likelihoods = [hom_ref_likelihood, het_likelihood, hom_alt_likelihood]
        gt = self.likelihoods_to_genotype(likelihoods)
        return SequenceCall.create(
            sequence=None,
            call_set=None,
            genotype=gt,
            genotype_likelihoods=likelihoods,
            info={
                "copy_number": float(
                    sequence_probe_coverage.median_depth) / expected_depth,
                "coverage": sequence_probe_coverage.coverage_dict,
                "expected_depths": self.expected_depths,
                "contamination_depths": self.contamination_depths})

    def _hom_alt_likeihood(self, median_depth, expected_depth):
        return log_lik_depth(median_depth, expected_depth * 0.75)

    def _het_likelihood(self, median_depth, expected_depth):
        return log_lik_depth(
            median_depth,
            expected_depth *
            self.minimum_detectable_frequency)

    def _hom_ref_likelihood(self, median_depth, expected_depth):
        return log_lik_depth(median_depth, expected_depth * 0.001)

    @property
    def minimum_detectable_frequency(self):
        if self.error_rate < 0.1:
            return 0.05
        else:
            return 0.25

    def _log_post_hom_ref(self, llk):
        return log(1) + llk

    def _log_post_het_or_alt(self, llk, expected_depth, sequence_coverage):
        expected_percentage_coverage = percent_coverage_from_expected_coverage(
            expected_depth)
        minimum_percentage_coverage_required = expected_percentage_coverage * \
            sequence_coverage.percent_coverage_threshold
        if sequence_coverage.percent_coverage > minimum_percentage_coverage_required:
            return self._log_post_hom_ref(llk)
        else:
            return MIN_LLK


class GeneCollectionTyper(Typer):

    """Types a collection of genes returning only the most likely version
        in the collection"""

    def __init__(self, depths, contamination_depths=[]):
        super(GeneCollectionTyper, self).__init__(depths, contamination_depths)
        self.presence_typer = PresenceTyper(depths, contamination_depths)

    def genotype(self, sequence_coverage_collection):
        """Types a collection of genes returning the most likely gene version
            in the collection with it's genotype"""
        best_version = self.get_best_version(
            sequence_coverage_collection.values())
        return self.presence_typer.genotype(best_version)

    def get_best_version(self, sequence_coverages):
        sequence_coverages.sort(key=lambda x: x.percent_coverage, reverse=True)
        current_best_gene = sequence_coverages[0]
        for gene in sequence_coverages[1:]:
            if gene.percent_coverage < current_best_gene.percent_coverage:
                return current_best_gene
            else:
                if gene.min_depth > current_best_gene.min_depth:
                    current_best_gene = gene
                elif gene.min_depth == current_best_gene.min_depth:
                    if gene.median_depth > current_best_gene.median_depth:
                        current_best_gene = gene
        return current_best_gene
