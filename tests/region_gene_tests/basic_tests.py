from atlas.gene import Region
from atlas.gene import Gene 
from nose.tools import assert_raises


class TestRegions():

    def setUp(self):
        pass

    def test_simple_gene(self):
        Gene(name = "gene", )
        
