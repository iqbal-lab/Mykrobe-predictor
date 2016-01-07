from unittest import TestCase
from atlas.typing import SequenceCoverage
from atlas.typing import PresenceTyper

class PresenceTyperTest(TestCase):

	def setUp(self):
		self.pt = PresenceTyper(depths = [100])
		self.pt_10 = PresenceTyper(depths = [10])

	def teardown(self):
		pass

	def test_base_case_no_coverage(self):
		s1 = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 0,
										median_depth = 0
										)
		vs = self.pt.genotype(s1)
		assert vs.gt == "0/0"

	def test_genotyping_gene_11(self):
		s = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 100,
										median_depth = 100,
										percent_coverage_threshold = 80,
										)
		s = self.pt.genotype(s)
		assert s.gt == "1/1"


	def test_genotyping_gene_01(self):
		s = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 82,
										median_depth = 2,
										percent_coverage_threshold = 80,
										)
		s = self.pt.genotype(s)
		assert s.gt == "0/1"	

	def test_resistotype_gene_at_high_CN(self):
		s = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 100,
										median_depth = 1000,
										percent_coverage_threshold = 80,
										)
		s = self.pt.genotype(s)
		assert s.gt == "1/1"

	def test_low_coverage(self):
		s = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 16,
										median_depth = 16,
										percent_coverage_threshold = 80,
										)
		s = self.pt_10.genotype(s)
		assert s.gt == "0/0"

		s = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 80,
										median_depth = 16,
										percent_coverage_threshold = 80,
										)
		s = self.pt_10.genotype(s)
		assert s.gt == "1/1"								

