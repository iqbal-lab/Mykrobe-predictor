from models import AlleleGenerator 
from models import Variant 
from nose.tools import assert_raises

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        pass

class TestSNPAlleleGenerator(BaseTest):

    def setUp(self):
        self.pg = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta")

    # def test_panel_generator(self):
    #     pg = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta")
    #     assert pg.ref is not None

    # def test_simple_variant(self):
    #     v = Variant("A", 31, "T")
    #     panel = self.pg.create(v)
    #     assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
    #     assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT"]

    #     v = Variant("A", 32, "T")
    #     panel = self.pg.create(v)
    #     assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
    #     assert panel.alts == ["GATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTGATC"]        

    # def test_simple_variant_invalid(self):
    #     with assert_raises(ValueError) as cm:
    #         v = Variant("T", 31, "T")
    #         panel = self.pg.create(v)

    # def test_simple_variant_start(self):
    #     v = Variant("C", 1, "T")
    #     panel = self.pg.create(v)
    #     assert panel.ref == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
    #     assert panel.alts == ["TGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"]

    # def test_simple_variant_end(self):
    #     v = Variant("A", 2902618, "T")     
    #     panel = self.pg.create(v)
    #     assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
    #     assert panel.alts == ["TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTTT"]

    #     v = Variant("T", 2902616, "C")
    #     panel = self.pg.create(v)
    #     assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
    #     assert panel.alts == ["TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTCTAT"]

    # def test_simple_variant_with_nearby_snp(self):
    #     v = Variant("A", 31, "T")
    #     v2 = Variant("A", 32, "T")
    #     panel = self.pg.create(v, context = [v2])
    #     assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
    #     assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
    #                           "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT"]    

    # def test_simple_variant_with_multiple_nearby_snps(self):
    #     v = Variant("A", 31, "T")
    #     v2 = Variant("A", 32, "T")
    #     v3 = Variant("C", 30, "G")

    #     panel = self.pg.create(v, context = [v2, v3])
    #     print panel.alts
    #     assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
    #     assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
    #                           "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT",
    #                           "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTGAT",
    #                           "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTGAT"]

    # def test_multiple_contexts(self):
    #     self.pg._create_multiple_contexts([v2, v3, v4, v5])

    def test_simple_variant_with_multiple_nearby_snps2(self):
        v = Variant("A", 31, "T")
        v2 = Variant("A", 32, "T")
        v3 = Variant("C", 30, "G")
        v4 = Variant("C", 30, "T")
        v5 = Variant("C", 30, "A")

        panel = self.pg.create(v, context = [v2, v3, v4, v5])
        print panel.alts
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert sorted(panel.alts) == sorted(["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGTTATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGTTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGATATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGATTTCAAATTTCATAACATCACCATGAGTTTGAT"])




