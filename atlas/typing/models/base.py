class ProbeCoverage(object):

    "Summary of kmer coverage of sequence. e.g output of color covearges"

    def __init__(self, percent_coverage, median_depth, min_depth):
        self.percent_coverage = percent_coverage
        self.median_depth = median_depth
        self.min_depth = min_depth