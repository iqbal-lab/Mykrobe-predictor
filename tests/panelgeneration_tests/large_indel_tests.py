from atlas.panelgeneration import AlleleGenerator 
from atlas.panelgeneration import Variant 
from nose.tools import assert_raises

class TestLargeINDELAlleleGenerator():

    def setUp(self):
        self.pg = AlleleGenerator(reference_filepath = "data/NC_000962.2.fasta")

    def test_large_variant(self):
        v = Variant("AACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACT", 1355983,"ACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCG")
        panel = self.pg.create(v)
        assert panel.ref ==   "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"
        assert panel.alts ==  ["TCGTCACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCGGGACC"]

    def test_large_variant2(self):
        v = Variant("AACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCAC", 1355983,"ACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGC")
        panel = self.pg.create(v)
        assert panel.ref ==   "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"
        assert panel.alts ==  ["TCGTCACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCTGGACC"] 


    def test_large_variant3(self):
        v = Variant("TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC", 1355978,"TCGTCAACGCCCGGTATCTGAGGATCGGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC")
        panel = self.pg.create(v)
        assert panel.ref ==   "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"
        assert panel.alts ==  ["TCGTCAACGCCCGGTATCTGAGGATCGGTGTTCTCACCCAATACAAGTCGCATTCACTGGACC"]                   