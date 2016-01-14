from unittest import TestCase
from atlas.typing import TypedVariant
from atlas.typing import VariantTyper

class VariantTyperTest(TestCase):

	def setUp(self):
		self.vt = VariantTyper(depths = [100])

	def teardown(self):
		pass

	def test_wt_vars(self):
		v1 = TypedVariant.create_object(name="A123T",
										reference_percent_coverage = 100,
										alternate_percent_coverage = 3,
										reference_median_depth = 100,
										alternate_median_depth = 100
										)
		vs = self.vt.type({"v1" : [v1]})
		assert vs.get("v1")[0].gt == "0/0"

	def test_alt_vars(self):
		v1 = TypedVariant.create_object(name="A123T",
										reference_percent_coverage = 3,
										alternate_percent_coverage = 100,
										reference_median_depth = 100,
										alternate_median_depth = 100
										)
		vs = self.vt.type({"v1" : [v1]})
		assert vs.get("v1")[0].gt == "1/1"


	def test_mixed_vars(self):
		v1 = TypedVariant.create_object(name="A123T",
										reference_percent_coverage = 100,
										alternate_percent_coverage = 100,
										reference_median_depth = 50,
										alternate_median_depth = 50
										)
		vs = self.vt.type({"v1" : [v1]})
		assert vs.get("v1")[0].gt == "0/1"	


class VariantTyperWithContamination(TestCase):

	def setUp(self):
		self.vt_no_contaim = VariantTyper(depths = [100], contamination_depths = [])
		self.vt_contaim = VariantTyper(depths = [80], contamination_depths = [20])

	def teardown(self):
		pass	

	def test_simple_case(self):
		v1 = TypedVariant.create_object(name="A123T",
										reference_percent_coverage = 100,
										alternate_percent_coverage = 100,
										reference_median_depth = 80,
										alternate_median_depth = 20
										)
		vs = self.vt_no_contaim.type({"v1" : [v1]})
		assert vs.get("v1")[0].gt == "0/1"

		vs = self.vt_contaim.type({"v1" : [v1]})
		assert vs.get("v1")[0].gt == "0/0"		

