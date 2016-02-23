import datetime
import json


class ProbeCoverage(object):

    def __init__(self, reference_percent_coverage,
                  alternate_percent_coverage,
                  reference_median_depth,
                  alternate_median_depth,
                  allele_name = None,
                  params = {}):
        if reference_median_depth is None:
            reference_median_depth = 0
        if alternate_median_depth is None:
            alternate_median_depth = 0
        self.reference_percent_coverage = reference_percent_coverage
        self.alternate_percent_coverage = alternate_percent_coverage
        self.reference_median_depth = reference_median_depth
        self.alternate_median_depth = alternate_median_depth
        self.allele_name = allele_name
        self.params = params

    @property
    def coverage_dict(self):
        return {"reference_percent_coverage" : self.reference_percent_coverage, 
            "alternate_percent_coverage" : self.alternate_percent_coverage,
            "reference_median_depth" : self.reference_median_depth,
            "alternate_median_depth" : self.alternate_median_depth}

