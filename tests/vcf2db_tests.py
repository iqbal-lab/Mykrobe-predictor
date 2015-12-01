from atlas.vcf2db import Variant 
from atlas.vcf2db import VariantSet
from atlas.vcf2db import Call 
from atlas.vcf2db import CallSet 
from atlas.vcf2db import Reference
from atlas.vcf2db import split_var_name
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
        cs = CallSet.create(sample_id = "C00123", name = "C00123")
        assert cs.name == "C00123"

class TestVariant(BaseTest):

    def setup(self):
        self.ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")
        self.vs = VariantSet.create(name = "C00123")

    def test_create_SNP(self):
        v1 = Variant.create(variant_set = self.vs, start = 0, end = 1,
                            reference_bases = "A", alternate_bases = ["T"],
                            reference = self.ref)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alt == "T"
        assert v1.length == 0


    def test_create_insertion(self):
        v1 = Variant.create(variant_set = self.vs,
                            start = 0, end = 1, reference_bases = "T",
                            alternate_bases = ["TA"], reference = self.ref)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alt == "TA"        
        assert v1.is_insertion == True
        assert v1.is_deletion == False
        assert v1.is_indel == True
        assert v1._length == None
        assert v1.length == 1
        assert v1._length == 1


    def test_create_deletion(self):
        v1 = Variant.create(variant_set = self.vs,
                            start = 0, end = 1, reference_bases = "AA",
                            alternate_bases = ["A"], reference = self.ref)
        assert v1.start == 0
        assert v1.end == 1
        assert v1.alt == "A"        
        assert v1.reference_bases == "AA"        
        assert v1.is_insertion == False        
        assert v1.is_deletion == True
        assert v1.is_indel == True
        assert v1.length == 1


    def test_split_name(self):
        name = "A12T"
        r,pos,a =  split_var_name(name)
        assert r == "A"
        assert pos == 12
        assert a == "T"

    def test_split_name_del(self):
        name = "AA12T"
        r,pos,a =  split_var_name(name)
        assert r == "AA"
        assert pos == 12
        assert a == "T"

    def test_split_name_ins(self):
        name = "A12TT"
        r,pos,a =  split_var_name(name)
        assert r == "A"
        assert pos == 12
        assert a == "TT" 

    def test_split_name2(self):
        name = "A12T/A"
        r,pos,a =  split_var_name(name)
        assert r == "A"
        assert pos == 12
        assert a == "T/A"   

    def test_split_name3(self):
        name = "C-54T"
        r,pos,a =  split_var_name(name)
        assert r == "C"
        assert pos == -54
        assert a == "T"                        

class TestCall(BaseTest):

    def setup(self):
        self.ref = Reference.create(name = "ref", length = 10000, source_accessions = "SRA_ABC123")
        self.vs = VariantSet.create(name = "C00123")  
        self.v1 = Variant.create(variant_set = self.vs, start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"], reference = self.ref)
        self.cs = CallSet.create(sample_id = "C00123", name = "C00123")  
        

    def test_create_SNP_call(self):
        c1 = Call.create(variant = self.v1, call_set = self.cs, genotype = "0/1", genotype_likelihood = 0.91)
        assert c1.call_set_name == "C00123"
        assert self.v1.call == c1
        assert c1.genotype == [0, 1]




      