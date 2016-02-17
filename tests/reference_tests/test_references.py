from atlas.references.models import Reference
from atlas.references.models import ReferenceSet
from mongoengine import connect
DB = connect('atlas-test')

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        DB.drop_database('atlas-test')


class TestReferenceSet(BaseTest):

    def test_create_reference_set(self):
    	reference =  ReferenceSet().create_and_save(name = "ref_set")
    	assert reference.name == "ref_set"
    	r = ReferenceSet.objects.get(name = "ref_set")
    	assert r == reference	

class TestReference(BaseTest):

    def setUp(self):
        self.reference_set =  ReferenceSet.create_and_save(name = "ref_set")

    def test_create_reference(self):
    	reference =  Reference().create_and_save(name = "ref", md5checksum = "sre", reference_sets = [self.reference_set])
    	assert reference.name == "ref"
    	r = Reference.objects.get(name = "ref")
    	assert r == reference
        assert r.reference_sets[0].name == "ref_set"
