import datetime
from atlas.vcf import VCF
from atlas.variants import VariantSet
from atlas.variants import VariantSetMetadata
from atlas.variants import Variant
from atlas.variants import CallSet
from atlas.variants import Call
from mongoengine import connect
from atlas.references import ReferenceSet
from atlas.references import Reference

DB = connect('atlas-test')


class BaseTest():

    def setUp(self):
        DB.drop_database('atlas-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.reference = Reference().create_and_save(
            name="NC_000962.3",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])

    def teardown(self):
        DB.drop_database('atlas-test')


class TestAddNewVariantSet(BaseTest):

    def test_add_new_vcf_variant_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        # We create a global variant set as well as one for the individual VCF
        assert VariantSet.objects().count() == 2
        vs = VariantSet.objects()[0]
        assert vs.name == "test.vcf"


class TestAddNewVariantSetMetaData(BaseTest):

    def test_add_new_vcf_variant_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantSetMetadata.objects().count() >= 2
        assert VariantSetMetadata.objects(key="KMER").count() == 2


class TestAddNewCallSet(BaseTest):

    def test_add_new_call_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        # Only one callset but the callset should belong to multiple variant
        # sets
        assert CallSet.objects().count() == 1
        assert CallSet.objects()[0].created_at <= datetime.datetime.now()
        assert len(CallSet.objects()[0].variant_sets) == 2


class TestVariantsAndCalls(BaseTest):

    def test_add_add_variants_and_calls(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert Call.objects().count() == 21
        assert Variant.objects().count() == 21


class TestAddSecondVCF(BaseTest):

    def setUp(self):
        DB.drop_database('atlas-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")
        self.reference = Reference().create_and_save(
            name="NC_000962.3",
            md5checksum="sre",
            reference_sets=[
                self.reference_set])
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()

    def test_add_second_vcf_variant_set(self):
        # This VCF only has one Variant which is not in the first VCF
        vcf = VCF(
            f="tests/vcf_tests/test2.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantSet.objects().count() == 3
        assert CallSet.objects().count() == 2
        assert Call.objects().count() == 42
        assert Variant.objects().count() == 22


class TestAddVCFwithIndels(BaseTest):

    def test_add_second_vcf_variant_set(self):
        # This VCF only has one Variant which is not in the first VCF
        vcf = VCF(
            f="tests/vcf_tests/test3.vcf",
            reference_set_id=self.reference_set.id,
            method="CORTEX")
        vcf.add_to_database()
        assert VariantSet.objects().count() == 2
        assert CallSet.objects().count() == 1
        assert Call.objects().count() == 106
        assert Variant.objects().count() == 106
        assert Variant.snps().count() == 89
        assert Variant.indels().count() == 17
        assert Variant.insertions().count() == 8
        assert Variant.deletions().count() == 8
        assert Variant.ph_snps.count() == 1
