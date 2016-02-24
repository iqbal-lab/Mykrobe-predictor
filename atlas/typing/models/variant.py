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
        return self.alternate_coverages[0]

    @property
    def coverage_dict(self):
        return {"reference_percent_coverage": self.reference_percent_coverage,
                "alternate_percent_coverage": self.alternate_percent_coverage,
                "reference_median_depth": self.reference_median_depth,
                "alternate_median_depth": self.alternate_median_depth,
                "reference_min_depth": self.reference_min_depth,
                "alternate_min_depth": self.alternate_min_depth,                
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

    @property 
    def alternate_median_depth(self):
        return self.best_alternate_coverage.median_depth   

    @property 
    def alternate_min_depth(self):
        return self.best_alternate_coverage.min_depth                       

