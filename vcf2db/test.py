from variants import Variant 
from variants import Call 
from mongoengine import connect
connect('atlas-test')

class TestVariant():

    def setup(self):
    	pass
        
    def tear_down(self):
    	pass

    def test_create_SNP(self):
        v1 = Variant.create(start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"])
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alt == "T"
        

class TestCall():

    def setup(self):
    	pass
        
    def tear_down(self):
    	pass

    def test_create_SNP_call(self):
        v1 = Variant.create(start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"])    	
        c1 = Call.create(variant = v1, call_set_name = "C00123", genotype = "0/1", genotype_likelihood = 0.91)
        assert c1.call_set_name == "C00123"
        assert v1.call == c1
        assert c1.genotype == [0, 1]

        