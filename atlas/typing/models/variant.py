import datetime
import json


class VariantProbeCoverage(object):

    def __init__(self, reference_coverage,
                 alternate_coverages,
                 var_name=None,
                 params={}):
        self.reference_coverage = reference_coverage
        self.alternate_coverages = alternate_coverages
        self.var_name = var_name
        self.params = params
        self.best_alternate_coverage = self._choose_best_alternate_coverage()

    def _choose_best_alternate_coverage(self):
        self.alternate_coverages.sort(
            key=lambda x: x.percent_coverage,
            reverse=True)
        current_best = self.alternate_coverages[0]
        for probe_coverage in self.alternate_coverages[1:]:
            if probe_coverage.percent_coverage < current_best.percent_coverage:
                return current_best
            else:
                if probe_coverage.min_depth > current_best.min_depth:
                    current_best = probe_coverage
                elif probe_coverage.min_depth <= current_best.min_depth:
                    if probe_coverage.median_depth > current_best.median_depth:
                        current_best = probe_coverage
        return current_best

    @property
    def coverage_dict(self):
        return {"reference": self.reference_coverage.coverage_dict,
                "alternate": self.best_alternate_coverage.coverage_dict
                }

    @property
    def reference_percent_coverage(self):
        return self.reference_coverage.percent_coverage

    @property
    def reference_median_depth(self):
        return self.reference_coverage.median_depth

    @property
    def reference_min_depth(self):
        return self.reference_coverage.min_depth

    @property
    def alternate_percent_coverage(self):
        return self.best_alternate_coverage.percent_coverage

    @alternate_percent_coverage.setter
    def alternate_percent_coverage(self, value):
        self.best_alternate_coverage.percent_coverage = value

    @property
    def alternate_median_depth(self):
        return self.best_alternate_coverage.median_depth

    @property
    def alternate_min_depth(self):
        return self.best_alternate_coverage.min_depth
