from atlas.panelgeneration import AlleleGenerator 
from atlas.panelgeneration import Variant 
from nose.tools import assert_raises
class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        pass

class TestSNPAlleleGenerator(BaseTest):

    def setUp(self):
        self.pg = AlleleGenerator(reference_filepath = "data/R00000022.fasta")
        # print self.pg.ref_length

    def test_panel_generator(self):
        pg = AlleleGenerator(reference_filepath = "data/R00000022.fasta")
        assert pg.ref is not None

    def test_simple_variant(self):
        v = Variant("A", 31, "T")
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT"]
        assert self.pg._calculate_length_delta_from_indels(v, []) == 0
        assert v.is_indel is False
        v = Variant("A", 32, "T")
        panel = self.pg.create(v)
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["GATTAAAGATAGAAATACACGATGCGAGCATTCAAATTTCATAACATCACCATGAGTTTGATC"]        

    def test_simple_variant_invalid(self):
        with assert_raises(ValueError) as cm:
            v = Variant("T", 31, "T")
            panel = self.pg.create(v)

    def test_simple_variant_start(self):
        v = Variant("C", 1, "T")
        panel = self.pg.create(v)
        assert panel.ref == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["TGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_variant_end(self):
        v = Variant("A", 2902618, "T")     
        panel = self.pg.create(v)
        assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert panel.alts == ["TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTTT"]

        v = Variant("T", 2902616, "C")
        panel = self.pg.create(v)
        assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert panel.alts == ["TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTCTAT"]

    def test_simple_variant_with_nearby_snp(self):
        v = Variant("A", 31, "T")
        v2 = Variant("A", 32, "T")
        panel = self.pg.create(v, context = [v2])
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT"]    

    def test_simple_variant_with_multiple_nearby_snps(self):
        v = Variant("A", 31, "T")
        v2 = Variant("A", 32, "T")
        v3 = Variant("C", 30, "G")

        panel = self.pg.create(v, context = [v2, v3])
        # print panel.alts
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTGAT",
                              "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_variant_with_multiple_nearby_snps2(self):
        v = Variant("A", 31, "T")
        v2 = Variant("A", 32, "T")
        v3 = Variant("C", 30, "G")
        v4 = Variant("C", 30, "T")
        v5 = Variant("C", 30, "A")
        panel = self.pg.create(v, context = [v2, v3, v4, v5])
        # for context in self.pg._create_multiple_contexts([v2, v3, v4, v5]):
        #     print [ "".join([v.ref, str(v.pos), v.alt]) for v in context]
        assert panel.ref ==                   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert sorted(panel.alts) == sorted([ "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGTTATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGTTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGATATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGATTTCAAATTTCATAACATCACCATGAGTTTGAT"])

    def test_simple_variant_with_multiple_nearby_snps(self):
        v = Variant("A", 31, "T")
        v2 = Variant("A", 32, "T")
        v5 = Variant("A", 32, "G")                        
        v3 = Variant("C", 30, "G")
        v4 = Variant("C", 30, "T")

        panel = self.pg.create(v, context = [v2, v3, v4, v5])
        # print sorted([v for v in self.pg._create_multiple_contexts([v2, v3, v4, v5])])
        # assert sorted([v for v in self.pg._create_multiple_contexts([v2, v3, v4, v5])]) == sorted([[v4, v5],
        #                                                                                           [v3, v5],
        #                                                                                           [v4, v2],
        #                                                                                           [v3, v2]])
        assert panel.ref ==                   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert sorted(panel.alts) == sorted([ "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGCTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGGTATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGGTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGTTATCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGTTTTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGCTGTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGGTGTCAAATTTCATAACATCACCATGAGTTTGAT",
                                              "CGATTAAAGATAGAAATACACGATGCGAGTTGTCAAATTTCATAACATCACCATGAGTTTGAT"
                                              ])

class TestINDELAlleleGenerator(TestSNPAlleleGenerator):


    def test_simple_deletion1(self):
        v = Variant("AA", 31, "A")
        assert v.is_indel
        assert v.is_deletion
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert self.pg._calculate_length_delta_from_indels(v, []) == 1
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTGATC"]

    def test_simple_deletion2(self):
        v = Variant("AT", 32, "A")
        panel = self.pg.create(v)
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATCC"
        assert panel.alts == ["GATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGATCC"]        

    def test_simple_deletion3(self):
        v = Variant("AT", 2902618, "T")
        panel = self.pg.create(v)
        assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "ATAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"        
        assert panel.alts == ["ATAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTT"] 

#         ## To do. 3 letter deletion

    def test_simple_insertion1(self):
        v = Variant("C", 1, "TTTC")
        panel = self.pg.create(v)
        assert v.is_indel
        assert v.is_insertion
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTT"        
        assert panel.alts == ["TTTCGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTT"]

    def test_simple_insertion2(self):
        v = Variant("C", 1, "CTTT")
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTT"                
        assert panel.alts == ["CTTTGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTT"]

    def test_simple_insertion3(self):
        v = Variant("A", 31, "ATTT")
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "ATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"                        
        assert panel.alts == ["ATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTGA"]

    def test_simple_insertion4(self):
        v = Variant("A", 32, "AGGGG")
        panel = self.pg.create(v)
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "TTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGA"                        
        assert panel.alts == ["TTAAAGATAGAAATACACGATGCGAGCAAGGGGTCAAATTTCATAACATCACCATGAGTTTGA"]        

    def test_simple_insertion5(self):     
        v = Variant("A", 2902618, "ATGC")
        panel = self.pg.create(v)
        assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert "".join(self.pg._get_alternate_reference_segment(v, [])) == "CAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"                
        assert panel.alts == ["CAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTATGCT"] 



