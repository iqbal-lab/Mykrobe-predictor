from models import Variant 
from models import VariantSet
from models import Call 
from models import CallSet 
from models import Reference
from mongoengine import connect
DB = connect('atlas-test')

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        DB.drop_database('atlas-test')

class TestReference(BaseTest):

    def test_create_SNP_call(self):
        ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")    	
        assert ref.name == "ref"

class TestVariantSet(BaseTest):

    def test_create_variant_set(self):
        vs = VariantSet.create(name = "C001234")        
        assert vs.name == "C001234"

class TestCallSet(BaseTest):

    def test_create_SNP_call(self):
        cs = CallSet.create(name = "C00123")    	
        assert cs.name == "C00123"

class TestVariant(BaseTest):

    def setup(self):
    	self.ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")    	    	
        self.vs = VariantSet.create(name = "C00123")  

    def test_create_SNP(self):
        v1 = Variant.create(variant_set = self.vs, start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"], reference = self.ref)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alt == "T"
        

class TestCall(BaseTest):

    def setup(self):
    	self.ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")    	    	
        self.vs = VariantSet.create(name = "C00123")          
        self.v1 = Variant.create(variant_set = self.vs, start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"], reference = self.ref)    	
        self.cs = CallSet.create(name = "C00123")  
        

    def test_create_SNP_call(self):
        c1 = Call.create(variant = self.v1, call_set = self.cs, genotype = "0/1", genotype_likelihood = 0.91)
        assert c1.call_set_name == "C00123"
        assert self.v1.call == c1
        assert c1.genotype == [0, 1]


        

      