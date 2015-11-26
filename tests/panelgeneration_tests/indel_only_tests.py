from atlas.panelgeneration import AlleleGenerator 
from atlas.panelgeneration import Variant 
from nose.tools import assert_raises


class TestINDELAlleleGenerator():

    def setUp(self):
        self.pg = AlleleGenerator(reference_filepath = "data/R00000022.fasta")
        self.pg2 = AlleleGenerator(reference_filepath = "data/NC_000962.2.fasta")

    def test_simple_deletion1(self):
        v = Variant("AA", 31, "A")
        assert v.is_indel
        assert v.is_deletion
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert self.pg._calculate_length_delta_from_indels(v, []) == 1
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTGATC"]

    def test_simple_deletion2(self):
        v = Variant("AT", 32, "A")
        panel = self.pg.create(v)
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["GATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGATCC"]        

    def test_simple_deletion3(self):
        v = Variant("AT", 2902618, "T")
        panel = self.pg.create(v)
        assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
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
        assert panel.alts == ["TTTCGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion2(self):
        v = Variant("C", 1, "CTTT")
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CTTTGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion3(self):
        v = Variant("A", 31, "ATTT")
        panel = self.pg.create(v)
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTGAT"]

    def test_simple_insertion4(self):
        v = Variant("A", 32, "AGGGG")
        panel = self.pg.create(v)
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["GATTAAAGATAGAAATACACGATGCGAGCAAGGGGTCAAATTTCATAACATCACCATGAGTTTGATC"]        

    def test_simple_insertion5(self):     
        v = Variant("A", 2902618, "ATGC")
        panel = self.pg.create(v)
        assert panel.ref ==   "TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTAT"
        assert panel.alts == ["TAACAAAATCCTTTTTATAACGCAAGTTCATTTTATACTACTGCTCAATTTTTTTACTTTTATGCT"] 

    def test_double_insertion(self):     
        v = Variant("A", 4021408, "ACGCTGGCGGGCG")
        v1  = Variant("AGA", 4021406, "CGG")
        context = [v1]
        assert self.pg2._remove_overlapping_contexts(v, [v1]) == []
        panel = self.pg2.create(v, context = context)
        assert panel.ref ==   "ATCTAGCCGCAAGGGCGCGAGCAGACGCAGAATCGCATGATTTGAGCTCAAATCATGCGATTC"
        assert panel.alts == ["ATCTAGCCGCAAGGGCGCGAGCAGACGCAGACGCTGGCGGGCGATCGCATGATTTGAGCTCAAATCATGCGATTC"] 

    def test_large_insertion(self):
        v = Variant("CCGCCGGCCCCGCCGTTT",1636155,"CTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCG")
        panel = self.pg2.create(v, context = [])
        assert panel.ref ==   "GAGTCGCCGAGGACGCCGGCGCCGCCATTGTCGCCAAATACCGTGAGACCTAGCAGGGTGCCGGCGCCGCCCTTGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCATCGCCGATGATGTTTTCCCCGCCCTTGCCGCCAGCCCCAGCGTTCCCG"        
        assert panel.alts == ["GGTTGGATCGCCACCGGCGCCACCGGCGCCGCCCGCGCCACCAGCACCGCCGCTGCCATCTGGGTCCGTCGAGTCGCCGAGGACGCCGGCGCCGCCATTGTCGCCAAATACCGTGAGACCTAGCAGGGTGCCGGCGCCGCCCTTGCTGCCGGCCCCGCCGGCGCCGCCCAATCCACCGAAGCCCCTCCCTTCGGTGGGGTCGCTGCCGCCGTCGCCGCCGTCACCGCCCTTGCCGCCGGCCCCGCCGTCGCCGCCGGCTCCGGCGGTGCCGTCGCCGCCCTGGCCGCCGGCCCCGCCGTTTCCGCCGCCGCCGCCATCGCCGATGATGTTTTCCCCGCCCTTGCCGCCAGCCCCAGCGTTCCCGCCGGCTCCGCCACTGGCGCCGGTGCCGCCGGGTGCAACGGCGTTGGCGCCGTTACCGCCGTTGCCGCCTTT"]





  # 