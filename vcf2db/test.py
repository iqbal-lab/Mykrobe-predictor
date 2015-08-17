from variants import Variant 
from variants import Call 
from variants import CallSet 
from references import Reference
from mongoengine import connect
connect('atlas-test')

class TestReference():

    def setup(self):
    	pass
        
    def tear_down(self):
    	pass

    def test_create_SNP_call(self):
        ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")    	
        assert ref.name == "ref"


class TestCallSet():

    def setup(self):
    	pass
        
    def tear_down(self):
    	pass

    def test_create_SNP_call(self):
        cs = CallSet.create(name = "C00123")    	
        assert cs.name == "C00123"

class TestVariant():

    def setup(self):
    	self.ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")    	    	

        
    def tear_down(self):
    	pass

    def test_create_SNP(self):
        v1 = Variant.create(start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"], reference = self.ref)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alt == "T"
        

class TestCall():

    def setup(self):
    	self.ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")    	    	
        self.v1 = Variant.create(start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"], reference = self.ref)    	
        self.cs = CallSet.create(name = "C00123")  
        
    def tear_down(self):
    	pass

    def test_create_SNP_call(self):
        c1 = Call.create(variant = self.v1, call_set = self.cs, genotype = "0/1", genotype_likelihood = 0.91)
        assert c1.call_set_name == "C00123"
        assert self.v1.call == c1
        assert c1.genotype == [0, 1]


        

      