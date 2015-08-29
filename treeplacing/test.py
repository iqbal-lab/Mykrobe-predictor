from unittest import TestCase

from models import Tree 
from models import Node 
from models import Leaf 
from models import Placer 
from nose.tools import assert_raises

from .. import Reference
from .. import VariantSet
from .. import Variant
from .. import CallSet
from .. import Call
from .. import GenotypedVariant

from mongoengine import connect
from pymongo import MongoClient
c = MongoClient()

DBNAME='atlas-test'
DB = connect(DBNAME)

class BaseTest():

    def setUp(self):
        pass
        
    def teardown(self):
        DB.drop_database(DBNAME)
        c.drop_database(DBNAME)
        Reference.drop_collection()



class TestNodes(BaseTest):

    def test_single_node_no_children(self):
        node = Node()
        assert node.children == []
        assert node.num_samples == 0
        assert node.is_leaf == False


    def test_node_triplet(self):
        node1 = Leaf(sample = 'C1')
        node2 = Leaf(sample = 'C2')
        root = Node(children = [node1, node2]) 
        assert root.num_samples == 2
        assert root.samples == ['C1', 'C2']

class TestMultiNode(TestNodes):

    def setUp(self):
        self.l1 = Leaf(sample = 'C1') 
        self.l2 = Leaf(sample = 'C2') 
        self.l3 = Leaf(sample = 'C3') 
        self.l4 = Leaf(sample = 'C4') 
        self.l5 = Leaf(sample = 'C5')         

        self.node1 = Node(children = [self.l1, self.l2])
        self.node2 = Node(children = [self.l4, self.l5])
        self.node3 = Node(children = [self.node2, self.l3])
        self.root = Node(children = [self.node1, self.node3])
        self.ref = Reference.create(name = "ref1", length = 10000, source_accessions = "SRA_ABC123")
        vs1 = VariantSet.create(name = "C1")
        cs1 = CallSet.create(name = "C1") 
        vs2 = VariantSet.create(name = "C2")
        cs2 = CallSet.create(name = "C2") 
        vs3 = VariantSet.create(name = "C3")
        cs3 = CallSet.create(name = "C3") 
        vs4 = VariantSet.create(name = "C4")
        cs4 = CallSet.create(name = "C4") 
        vs5 = VariantSet.create(name = "C5")
        cs5 = CallSet.create(name = "C5")                 


        self.v1 = Variant.create(variant_set = vs1,
                                 start = 1,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v2 = Variant.create(variant_set = vs2,
                                 start = 2,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v3 = Variant.create(variant_set = vs3,
                                 start = 3,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v4 = Variant.create(variant_set = vs4,
                                 start = 4,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v5 = Variant.create(variant_set = vs5,
                                 start = 4,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)                                                                                                           

        Call.create(variant = self.v1, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v2, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v3, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v4, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v5, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)


    def test_multi_node(self):
        assert self.root.num_samples == 5
        assert sorted(self.root.samples) == ['C1', 'C2', 'C3', 'C4', 'C5']

    def test_phylo_snps(self):
        assert self.node1.phylo_snps == [self.v1.name, self.v2.name]
        assert self.node2.phylo_snps == [self.v4.name]
        assert self.node3.phylo_snps == [self.v3.name, self.v4.name]
        assert self.root.phylo_snps == [self.v1.name, self.v2.name, self.v3.name, self.v4.name]

        assert self.l1.phylo_snps == [self.v1.name]        
        assert self.l2.phylo_snps == [self.v2.name]        
        assert self.l3.phylo_snps == [self.v3.name]
        assert self.l4.phylo_snps == []        
        assert self.l5.phylo_snps == []

    def test_placement(self):
        new_call_set = CallSet.create(name = "C6") 
        GenotypedVariant.create("A1T", new_call_set, 30)
        assert Placer(root = self.root).place("C6") == "C1"

    def test_abigious_placement(self):
        new_call_set = CallSet.create(name = "C7") 
        GenotypedVariant.create("A4T", new_call_set, 30)
        assert Placer(root = self.root).place("C7") == ["C4", "C5"]


class TestMultiNodeHomoplasy(TestNodes):

    def setUp(self):
        self.l1 = Leaf(sample = 'C1') 
        self.l2 = Leaf(sample = 'C2') 
        self.l3 = Leaf(sample = 'C3') 
        self.l4 = Leaf(sample = 'C4') 
        self.l5 = Leaf(sample = 'C5')         

        self.node1 = Node(children = [self.l1, self.l2])
        self.node2 = Node(children = [self.l4, self.l5])
        self.node3 = Node(children = [self.node2, self.l3])
        self.root = Node(children = [self.node1, self.node3])
        self.ref = Reference.create(name = "ref1", length = 10000, source_accessions = "SRA_ABC123")
        vs1 = VariantSet.create(name = "C1")
        cs1 = CallSet.create(name = "C1") 
        vs2 = VariantSet.create(name = "C2")
        cs2 = CallSet.create(name = "C2") 
        vs3 = VariantSet.create(name = "C3")
        cs3 = CallSet.create(name = "C3") 
        vs4 = VariantSet.create(name = "C4")
        cs4 = CallSet.create(name = "C4") 
        vs5 = VariantSet.create(name = "C5")
        cs5 = CallSet.create(name = "C5")                 


        self.v1 = Variant.create(variant_set = vs1,
                                 start = 1,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v2 = Variant.create(variant_set = vs2,
                                 start = 2,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v3 = Variant.create(variant_set = vs3,
                                 start = 3,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v4 = Variant.create(variant_set = vs4,
                                 start = 4,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)

        self.v5 = Variant.create(variant_set = vs5,
                                 start = 1,
                                 end = 2,
                                 reference_bases = "A",
                                 alternate_bases = ["T"],
                                 reference = self.ref)                                                                                                           

        Call.create(variant = self.v1, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v2, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v3, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v4, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)
        Call.create(variant = self.v5, call_set = cs1, genotype = "0/1", genotype_likelihood = 0.91)


    def test_multi_node(self):
        assert self.root.num_samples == 5
        assert sorted(self.root.samples) == ['C1', 'C2', 'C3', 'C4', 'C5']

    def test_phylo_snps(self):
        assert self.node1.phylo_snps == [self.v2.name]
        assert self.node2.phylo_snps == [self.v4.name]
        assert self.node3.phylo_snps == [self.v3.name, self.v4.name]
        assert self.root.phylo_snps == [self.v1.name, self.v2.name, self.v3.name, self.v4.name]

        assert self.l1.phylo_snps == [self.v1.name]        
        assert self.l2.phylo_snps == [self.v2.name]        
        assert self.l3.phylo_snps == [self.v3.name]
        assert self.l4.phylo_snps == []        
        assert self.l5.phylo_snps == []

    def test_placement(self):
        new_call_set = CallSet.create(name = "C6") 
        GenotypedVariant.create("A1T", new_call_set, 30)
        assert Placer(root = self.root).place("C6") == "C1"

    def test_abigious_placement(self):
        new_call_set = CallSet.create(name = "C7") 
        GenotypedVariant.create("A4T", new_call_set, 30)
        assert Placer(root = self.root).place("C7") == ["C4", "C5"]




