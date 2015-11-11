from atlas.panelgeneration import AlleleGenerator 
from atlas.panelgeneration import Variant 
from nose.tools import assert_raises


class TestINDELandSNPSAlleleGenerator():

    def setUp(self):
        self.pg = AlleleGenerator(reference_filepath = "data/R00000022.fasta")

    def test_ins_with_SNP_context(self):
        v = Variant("A", 31, "ATTT")
        v2 = Variant("A", 32, "T")
        panel = self.pg.create(v, context = [v2])
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["ATTAAAGATAGAAATACACGATGCGAGCATTTATCAAATTTCATAACATCACCATGAGTTTGA",
                              "ATTAAAGATAGAAATACACGATGCGAGCATTTTTCAAATTTCATAACATCACCATGAGTTTGA"]

    def test_del_with_SNP_context1(self):
        v = Variant("AA", 31, "A")
        v2 = Variant("T", 33, "A")
        panel = self.pg.create(v, context = [v2])
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTGATC",
        					  "CGATTAAAGATAGAAATACACGATGCGAGCAACAAATTTCATAACATCACCATGAGTTTGATC"]                                 

    def test_del_with_SNP_context2(self):
        v = Variant("AA", 31, "A")
        v2 = Variant("A", 32, "T")
        panel = self.pg.create(v, context = [v2])
        self.pg._remove_contexts_spanning_del(v, [v2]) == []
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCATCAAATTTCATAACATCACCATGAGTTTGATC"]

    def test_del_with_ins_context1(self):
        v = Variant("AAT", 31, "A")
        v2 = Variant("T", 4, "TTTT")
        panel = self.pg.create(v, context = [v2])
        self.pg._remove_contexts_spanning_del(v, [v2]) == [v2]
        assert panel.ref ==   "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCACAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTTTTAAAGATAGAAATACACGATGCGAGCACAAATTTCATAACATCACCATGAGTTTGAT"]  


    def test_del_with_ins_context2(self):
        v = Variant("ATC", 32, "A")
        v2 = Variant("C", 1, "CTTT")
        panel = self.pg.create(v, context = [v2])
        assert self.pg._remove_contexts_spanning_del(v, [v2]) == [v2]
        assert self.pg._remove_contexts_not_within_k(v, [v2]) == []
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"] 

    def test_del_with_ins_context3(self):
        v = Variant("ATC", 32, "A")
        v2 = Variant("T", 5, "TT")
        panel = self.pg.create(v, context = [v2])
        assert self.pg._remove_contexts_spanning_del(v, [v2]) == [v2]
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"]


    def test_del_with_ins_context4(self):
        v = Variant("ATC", 32, "A")
        v2 = Variant("T", 5, "TT")
        v3 = Variant("T", 5, "TG")
        panel = self.pg.create(v, context = [v2, v3])
        assert self.pg._remove_contexts_spanning_del(v, [v2, v3]) == [v2, v3]
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert sorted(panel.alts) == sorted(["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTGAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"])



    def test_del_with_ins_context5(self):
        v = Variant("ATC", 32, "A")
        v2 = Variant("T", 5, "TT")
        v3 = Variant("A", 6, "AG")
        panel = self.pg.create(v, context = [v2, v3])
        assert self.pg._remove_contexts_spanning_del(v, [v2, v3]) == [v2, v3]
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTAGAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                              "GATTTAGAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATC"]


    def test_del_with_ins_context_where_base_is_deleted(self):
        v = Variant("ATC", 32, "A")
        v2 = Variant("T", 33, "C")
        panel = self.pg.create(v, context = [v2])
        assert self.pg._remove_contexts_spanning_del(v, [v2]) == []
        assert panel.ref ==   "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert panel.alts == ["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"]                                                                                                                

    def test_del_with_ins_context_where_base_is_deleted2(self):
        v = Variant("ATC", 32, "A")
        v2 = Variant("TAAA", 5, "T")
        v3 = Variant("A", 7, "AG")
        panel = self.pg.create(v, context = [v2, v3])
        assert panel.ref ==                  "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert sorted(panel.alts) == sorted(["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                                             "CGATTGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCCAAA",
                                             "GATTAAGAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"])

        panel = self.pg.create(v, context = [v3, v2])
        assert panel.ref ==                  "GATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATC"
        assert sorted(panel.alts) == sorted(["CGATTAAAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC",
                                             "CGATTGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCCAAA",
                                             "GATTAAGAGATAGAAATACACGATGCGAGCAAAAATTTCATAACATCACCATGAGTTTGATCC"])        

   	
