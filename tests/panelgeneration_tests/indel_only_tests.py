from atlas.panelgeneration import AlleleGenerator 
from atlas.panelgeneration import Variant 
from nose.tools import assert_raises


class TestINDELAlleleGenerator():

    def setUp(self):
        self.pg = AlleleGenerator(reference_filepath = "data/R00000022.fasta")
        # print self.pg.ref_length

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

    def test_simple_deletion4(self):
        v = Variant("ATC", 32, "A")
        panel = self.pg.create(v)
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"] 

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




  