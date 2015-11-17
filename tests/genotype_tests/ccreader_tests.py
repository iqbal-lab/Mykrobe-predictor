from atlas.genotyping import ColourCovgsReader
from nose.tools import assert_raises


class TestColourCovgsReader():

    def setUp(self):
    	self.f = open("tests/genotype_tests/test.fasta.colour_covgs", "r")

    def tearDown(self):
    	self.f.close()

    def test_basic_reader(self):

    	ccreader = ColourCovgsReader(self.f)

    	ccread = next(ccreader)
    	assert ccread.name == "ref-A1004177C"
    	assert ccread.median_non_zero_coverage == 0
    	assert ccread.percent_non_zero_coverage == 0

    	ccread2 = next(ccreader)
    	assert ccread2.name == "alt-A1004177C"
    	assert ccread2.percent_non_zero_coverage == 100 
    	assert ccread2.median_non_zero_coverage == 56

    	for i,read in enumerate(ccreader):
    		if i == 0:
    			assert read.name == "ref-A1024346G"
    		elif i == 1:
    			assert read.name == "alt-A1024346G"
    		elif i == 2:
    			assert read.name == "ref-A103600G"
    		elif i == 3:
    			assert read.name == "alt-A103600G"    			




