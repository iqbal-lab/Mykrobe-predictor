from atlas.variants import Variant 
from atlas.variants import VariantSet
from atlas.variants import Call 
from atlas.variants import CallSet 
from atlas.variants import Reference
from atlas.variants import split_var_name
from mongoengine import connect
DB = connect('atlas-test')

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        DB.drop_database('atlas-test')
