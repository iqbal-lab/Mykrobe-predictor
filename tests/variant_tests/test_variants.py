from atlas.variants import Variant 
from atlas.variants import VariantSet
from atlas.variants import Call 
from atlas.variants import CallSet 
from atlas.references import Reference
from atlas.references import ReferenceSet
from mongoengine import connect
DB = connect('atlas-test')

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        DB.drop_database('atlas-test')


class TestHighLevelDocuments(BaseTest):

    def test_create_reference_set(self):
    	reference =  ReferenceSet().create_and_save(name = "ref", md5checksum = "sre")
    	assert reference.name == "ref"
    	r = Reference.objects.get(name = "ref")
    	assert r == reference	

    def test_create_reference(self):
    	reference =  Reference().create_and_save(name = "ref", md5checksum = "sre")
    	assert reference.name == "ref"
    	r = Reference.objects.get(name = "ref")
    	assert r == reference

    def test_create_variant_set(self):
    	variant_set = VariantSet.create_and_save()

# class TestVariants(BaseTest):

#     def setUp(self):
#         reference =  Reference.create_and_save()

#     def test_create_new_variant(self):
#     	variant = Variant.create()
