from models import AlleleGenerator 
from nose.tools import assert_raises

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
    	pass

class TestReference(BaseTest):

    def setUp(self):
    	self.pg = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta")


    def test_panel_generator(self):
    	pg = AlleleGenerator(reference_filepath = "/home/phelimb/git/atlas/data/R00000022.fasta")
    	assert pg.ref is not None
    	assert len(pg.ref) > 100

    def test_simple_variant(self):
    	panel = self.pg.create("A", 31, "T")
    	assert str(panel.ref) == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
    	assert str(panel.alt) == "CGATTAAAGATAGAAATACACGATGCGAGCTATCAAATTTCATAACATCACCATGAGTTTGAT"

    def test_simple_variant_invalid(self):
    	with assert_raises(ValueError) as cm:
    		panel = self.pg.create("T", 31, "T")

    def test_simple_variant_start(self):
    	panel = self.pg.create("C", 1, "T")
    	assert str(panel.ref) == "CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"
    	assert str(panel.alt) == "TGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGAT"    		

