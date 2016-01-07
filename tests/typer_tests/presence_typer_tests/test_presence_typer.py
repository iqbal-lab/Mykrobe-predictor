from unittest import TestCase
from atlas.typing import SequenceCoverage
from atlas.typing import PresenceTyper

class PresenceTyperTest(TestCase):

	def setUp(self):
		self.pt = PresenceTyper(depths = [100])

	def teardown(self):
		pass

	def test_base_case_no_coverage(self):
		s1 = SequenceCoverage.create_object(name="A123T",
										percent_coverage = 0,
										median_depth = 0
										)
		vs = self.pt.genotype(s1)
		assert vs.gt == "0/0"

	# def test_alt_vars(self):
	# 	v1 = SequenceCoverage.create_object(name="A123T",
	# 									percent_coverage = 100,
	# 									median_depth = 100
	# 									)
	# 	vs = self.pt.type({"v1" : [v1]})
	# 	assert vs.get("v1")[0].gt == "1/1"


	# def test_mixed_vars(self):
	# 	v1 = SequenceCoverage.create_object(name="A123T",
	# 									percent_coverage = 100,
	# 									median_depth = 50
	# 									)
	# 	vs = self.pt.type({"v1" : [v1]})
	# 	assert vs.get("v1")[0].gt == "0/1"					
