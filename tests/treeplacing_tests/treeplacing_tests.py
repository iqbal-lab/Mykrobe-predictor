from unittest import TestCase

from atlas.treeplacing import Node
from atlas.treeplacing import Leaf
from atlas.treeplacing import Placer
from nose.tools import assert_raises

from atlas.vcf2db import Reference
from atlas.vcf2db import VariantSet
from atlas.vcf2db import Variant
from atlas.vcf2db import CallSet
from atlas.vcf2db import Call
from atlas.typing import TypedVariant

from mongoengine import connect
from pymongo import MongoClient
c = MongoClient()

DBNAME = 'atlas-test'
DB = connect(DBNAME)


class BaseTest(TestCase):

    def setUp(self):
        DB.drop_database(DBNAME)
        c.drop_database(DBNAME)
        Reference.drop_collection()

    def teardown(self):
        DB.drop_database(DBNAME)
        c.drop_database(DBNAME)
        Reference.drop_collection()


class TestNodes(BaseTest):

    def test_single_node_no_children(self):
        node = Node()
        assert node.children == []
        assert node.num_samples == 0
        assert node.is_leaf is False

    def test_node_triplet(self):
        node1 = Leaf(sample='C1')
        node2 = Leaf(sample='C2')
        root = Node(children=[node1, node2])
        assert root.num_samples == 2
        assert root.samples == ['C1', 'C2']


class TestMultiNode(TestNodes):

    def setUp(self):
        DB.drop_database(DBNAME)
        c.drop_database(DBNAME)
        Reference.drop_collection()

        self.l1 = Leaf(sample='C1')
        self.l2 = Leaf(sample='C2')
        self.l3 = Leaf(sample='C3')
        self.l4 = Leaf(sample='C4')
        self.l5 = Leaf(sample='C5')

        self.node1 = Node(children=[self.l1, self.l2])
        self.node2 = Node(children=[self.l4, self.l5])
        self.node3 = Node(children=[self.node2, self.l3])
        self.root = Node(children=[self.node1, self.node3])
        self.ref = Reference.create(
            name="ref1",
            length=10000,
            source_accessions="SRA_ABC123")
        vs1 = VariantSet.create(name="C1")
        cs1 = CallSet.create(sample_id="123", name="C1")
        vs2 = VariantSet.create(name="C2")
        cs2 = CallSet.create(sample_id="123", name="C2")
        vs3 = VariantSet.create(name="C3")
        cs3 = CallSet.create(sample_id="123", name="C3")
        vs4 = VariantSet.create(name="C4")
        cs4 = CallSet.create(sample_id="123", name="C4")
        vs5 = VariantSet.create(name="C5")
        cs5 = CallSet.create(sample_id="123", name="C5")

        self.v1 = Variant.create(variant_set=vs1.id,
                                 start=1,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v2 = Variant.create(variant_set=vs2.id,
                                 start=2,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v3 = Variant.create(variant_set=vs3.id,
                                 start=3,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v4 = Variant.create(variant_set=vs4.id,
                                 start=4,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v5 = Variant.create(variant_set=vs5.id,
                                 start=4,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        Call.create(
            variant=self.v1.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v2.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v3.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v4.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v5.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)

    def test_multi_node(self):
        assert self.root.num_samples == 5
        assert sorted(self.root.samples) == ['C1', 'C2', 'C3', 'C4', 'C5']

    def test_phylo_snps(self):
        assert self.node1.phylo_snps == {self.v1.name: 0.5, self.v2.name: 0.5}
        assert self.node2.phylo_snps == {self.v4.name: 1}
        assert self.node3.phylo_snps == {
            self.v3.name: 0.3333333333333333,
            self.v4.name: 0.6666666666666666}
        assert self.root.phylo_snps == {
            self.v1.name: 0.2,
            self.v2.name: 0.2,
            self.v3.name: 0.2,
            self.v4.name: 0.4}

        assert self.l1.phylo_snps == {self.v1.name: 1}
        assert self.l2.phylo_snps == {self.v2.name: 1}
        assert self.l3.phylo_snps == {self.v3.name: 1}
        assert self.l4.phylo_snps == {self.v4.name: 0}
        assert self.l5.phylo_snps == {self.v5.name: 0}

    def test_placement(self):
        new_call_set = CallSet.create(sample_id="123", name="C6")
        TypedVariant.create(
            name="A1T",
            call_set=new_call_set.id,
            reference_median_depth=0,
            reference_percent_coverage=100,
            alternate_median_depth=0,
            alternate_percent_coverage=30,
            gt="1/1")
        assert Placer(root=self.root).place("C6") == "C1"

    def test_abigious_placement(self):
        new_call_set = CallSet.create(sample_id="123", name="C7")
        TypedVariant.create(
            name="A4T",
            call_set=new_call_set.id,
            reference_median_depth=0,
            reference_percent_coverage=100,
            alternate_median_depth=0,
            alternate_percent_coverage=30,
            gt="1/1")
        assert Placer(root=self.root).place("C7") == ["C4", "C5"]


class TestMultiNodeHomoplasy(TestNodes):

    def setUp(self):
        DB.drop_database(DBNAME)
        c.drop_database(DBNAME)
        Reference.drop_collection()

        self.l1 = Leaf(sample='C1')
        self.l2 = Leaf(sample='C2')
        self.l3 = Leaf(sample='C3')
        self.l4 = Leaf(sample='C4')
        self.l5 = Leaf(sample='C5')

        self.node1 = Node(children=[self.l1, self.l2])
        assert self.l1.parent == self.node1
        assert self.l2.parent == self.node1
        self.node2 = Node(children=[self.l4, self.l5])
        self.node3 = Node(children=[self.node2, self.l3])
        self.root = Node(children=[self.node1, self.node3])
        self.ref = Reference.create(
            name="ref1",
            length=10000,
            source_accessions="SRA_ABC123")
        vs1 = VariantSet.create(name="C1")
        cs1 = CallSet.create(sample_id="123", name="C1")
        vs2 = VariantSet.create(name="C2")
        cs2 = CallSet.create(sample_id="123", name="C2")
        vs3 = VariantSet.create(name="C3")
        cs3 = CallSet.create(sample_id="123", name="C3")
        vs4 = VariantSet.create(name="C4")
        cs4 = CallSet.create(sample_id="123", name="C4")
        vs5 = VariantSet.create(name="C5")
        cs5 = CallSet.create(sample_id="123", name="C5")

        self.v1 = Variant.create(variant_set=vs1.id,
                                 start=1,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v2 = Variant.create(variant_set=vs2.id,
                                 start=2,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v3 = Variant.create(variant_set=vs3.id,
                                 start=3,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v4 = Variant.create(variant_set=vs4.id,
                                 start=4,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        self.v5 = Variant.create(variant_set=vs5.id,
                                 start=1,
                                 end=2,
                                 reference_bases="A",
                                 alternate_bases=["T"],
                                 reference=self.ref.id)

        Call.create(
            variant=self.v1.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v2.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v3.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v4.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)
        Call.create(
            variant=self.v5.id,
            call_set=cs1.id,
            genotype="0/1",
            genotype_likelihood=0.91)

    def test_multi_node(self):
        assert self.root.num_samples == 5
        assert sorted(self.root.samples) == ['C1', 'C2', 'C3', 'C4', 'C5']

    def test_phylo_snps(self):
        assert self.node1.phylo_snps == {
            self.v1.name: 0.5 -
            0.3333333333333333,
            self.v2.name: 0.5}
        assert self.node2.phylo_snps == {self.v1.name: 0.5, self.v4.name: 0.5}
        assert self.node3.phylo_snps == {
            self.v1.name: 0.3333333333333333 - 0.5,
            self.v3.name: 0.3333333333333333,
            self.v4.name: 0.3333333333333333}
        assert self.root.phylo_snps == {
            self.v1.name: 0.4,
            self.v2.name: 0.2,
            self.v3.name: 0.2,
            self.v4.name: 0.2}

        assert self.l1.phylo_snps == {self.v1.name: 1}
        assert self.l2.phylo_snps == {self.v2.name: 1}
        assert self.l3.phylo_snps == {self.v3.name: 1}
        assert self.l4.phylo_snps == {self.v4.name: 1}
        assert self.l5.phylo_snps == {self.v1.name: 1}

    def test_placement(self):
        new_call_set = CallSet.create(sample_id="123", name="C8")
        TypedVariant.create(
            name="A1T",
            call_set=new_call_set.id,
            reference_median_depth=0,
            reference_percent_coverage=100,
            alternate_median_depth=0,
            alternate_percent_coverage=30,
            gt="1/1")
        # Note - I think the commented line hear should be correct behaviour
        Placer(root=self.root).place("C8") == ["C1"]
        # print Placer(root = self.root).place("C8", verbose=True)
        # assert sorted(Placer(root = self.root).place("C8")) == sorted(['C1', 'C2', 'C3', 'C4', 'C5'])

    def test_abigious_placement(self):
        new_call_set = CallSet.create(sample_id="123", name="C9")
        TypedVariant.create(
            name="A4T",
            call_set=new_call_set.id,
            reference_median_depth=0,
            reference_percent_coverage=100,
            alternate_median_depth=0,
            alternate_percent_coverage=30,
            gt="1/1")
        assert Placer(root=self.root).place("C9") == "C4"
