from variants import Variant 
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
        
