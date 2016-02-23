from atlas.typing.typer.base import Typer
from atlas.stats import log_lik_R_S_coverage
from atlas.typing.typer.base import MIN_LLK

DEFAULT_ERROR_RATE = 0.05
DEFAULT_MINOR_FREQ = 0.1


class VariantTyper(Typer):

    def __init__(self, expected_depths, contamination_depths=[],
                 error_rate=DEFAULT_ERROR_RATE,
                 minor_freq=DEFAULT_MINOR_FREQ):
        super(
            VariantTyper,
            self).__init__(
            expected_depths,
            contamination_depths,
            error_rate)
        self.method = "MAP"
        self.error_rate = error_rate
        self.minor_freq = minor_freq

        if len(expected_depths) > 1:
            raise NotImplementedError("Mixed samples not handled yet")

    def type(self, probe_coverage):
        hom_ref_likelihood = self._hom_ref_lik(probe_coverage)
        hom_alt_likelihood = self._hom_alt_lik(probe_coverage)
        if not self.has_contamination():
            het_likelihood = self._het_lik(probe_coverage)
        else:
            het_likelihood = MIN_LLK
        gt = self.likelihoods_to_genotype(
            [hom_ref_likelihood, het_likelihood, hom_alt_likelihood])
        return {
            probe_coverage.allele_name: {
                "gt": gt,
                "coverage": probe_coverage.coverage_dict,
                "copy_number": float(
                    probe_coverage.alternate_median_depth) /
                self.expected_depths[0]}}

    def _hom_ref_lik(self, variant):
        if variant.reference_percent_coverage < 100:
            return MIN_LLK
        else:
            hom_ref_likes = []
            # Either alt+cov or alt_covg + contam_covg
            for expected_depth in self.expected_depths:
                hom_ref_likes.append(
                    log_lik_R_S_coverage(
                        variant.reference_median_depth,
                        variant.alternate_median_depth,
                        expected_depth,
                        expected_depth *
                        self.error_rate /
                        3))
                for contamination in self.contamination_depths:
                    hom_ref_likes.append(
                        log_lik_R_S_coverage(
                            variant.reference_median_depth,
                            variant.alternate_median_depth,
                            expected_depth + contamination,
                            (expected_depth + contamination) * self.error_rate / 3))
            return max(hom_ref_likes)

    def _hom_alt_lik(self, variant):
        if variant.alternate_percent_coverage < 100:
            return MIN_LLK
        else:
            hom_alt_liks = []
            # Either alt+cov or alt_covg + contam_covg
            for expected_depth in self.expected_depths:
                hom_alt_liks.append(
                    log_lik_R_S_coverage(
                        variant.alternate_median_depth,
                        variant.reference_median_depth,
                        expected_depth,
                        expected_depth *
                        self.error_rate /
                        3))
                for contamination in self.contamination_depths:
                    hom_alt_liks.append(
                        log_lik_R_S_coverage(
                            variant.alternate_median_depth,
                            variant.reference_median_depth,
                            expected_depth + contamination,
                            (expected_depth + contamination) * self.error_rate / 3))
            return max(hom_alt_liks)

    def _het_lik(self, variant):
        if variant.alternate_percent_coverage < 100 or variant.reference_percent_coverage < 100:
            return MIN_LLK
        else:
            het_liks = []
            for expected_depth in self.expected_depths:
                het_liks.append(
                    log_lik_R_S_coverage(
                        variant.alternate_median_depth,
                        variant.reference_median_depth,
                        expected_depth * self.minor_freq,
                        expected_depth * (
                            1 - self.minor_freq)))
            return max(het_liks)
