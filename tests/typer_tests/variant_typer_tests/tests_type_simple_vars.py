from unittest import TestCase
from atlas.schema import Variant
from atlas.schema import VariantCall
from atlas.typing import VariantTyper
from atlas.typing import ProbeCoverage
from atlas.typing import SequenceProbeCoverage
from atlas.typing import VariantProbeCoverage


class VariantTyperTest(TestCase):

    def setUp(self):
        self.vt = VariantTyper(expected_depths=[100])

    def teardown(self):
        pass

    def test_wt_vars(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=100)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=3,
                                             median_depth=100)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type([v1])
        assert call.genotype == [0, 0]
        assert call.info.get('expected_depths') == [100]

    def test_alt_vars(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=3,
                                           median_depth=100)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=100,
                                             median_depth=100)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type([v1])
        assert call.genotype == [1, 1]

    def test_mixed_vars(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=50)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=100,
                                             median_depth=50)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type(v1)
        assert call.genotype == [0, 1]


class VariantTyperWithContamination(TestCase):

    def setUp(self):
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[])
        self.vt_contaim = VariantTyper(
            expected_depths=[80],
            contamination_depths=[20])

    def teardown(self):
        pass

    def test_simple_case(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=80)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=100,
                                             median_depth=20)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )

        call = self.vt_no_contaim.type(v1)
        assert call.genotype == [0, 1]

        call = self.vt_contaim.type(v1)
        assert call.genotype == [0, 0]


class VariantTyperWithMultipleAlternateCoverages(TestCase):

    def setUp(self):
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[])

    def teardown(self):
        pass

    def test_simple_case(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=70,
                                           median_depth=80)
        alt1 = ProbeCoverage(min_depth=100,
                             percent_coverage=70,
                             median_depth=20)
        alt2 = ProbeCoverage(min_depth=100,
                             percent_coverage=100,
                             median_depth=80)
        alternate_coverages = [alt1, alt2]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )
        assert v1._choose_best_alternate_coverage() == alt2

        call = self.vt_no_contaim.type(v1)
        assert call.genotype == [1, 1]


class VariantTyperWithMultipleProbeCoverages(TestCase):

    def setUp(self):
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[])

    def teardown(self):
        pass

    def test_simple_case(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=80)
        alt1 = ProbeCoverage(min_depth=100,
                             percent_coverage=50,
                             median_depth=20)
        alt2 = ProbeCoverage(min_depth=100,
                             percent_coverage=40,
                             median_depth=80)
        alternate_coverages = [alt1, alt2]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )

        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=80,
                                           median_depth=80)
        alt1 = ProbeCoverage(min_depth=100,
                             percent_coverage=50,
                             median_depth=20)
        alt2 = ProbeCoverage(min_depth=100,
                             percent_coverage=100,
                             median_depth=80)

        alternate_coverages = [alt1, alt2]

        v2 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverage=reference_coverage,
                                  alternate_coverages=alternate_coverages
                                  )

        call = self.vt_no_contaim.type([v1, v2])
        assert call.genotype == [1, 1]     
