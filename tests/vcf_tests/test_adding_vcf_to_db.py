from atlas.vcf import VCF
from atlas.variants import VariantSet
from atlas.variants import VariantSetMetadata
from mongoengine import connect
from atlas.references import ReferenceSet

DB = connect('atlas-test')


class BaseTest():

    def setUp(self):
        DB.drop_database('atlas-test')
        self.reference_set = ReferenceSet().create_and_save(name="ref_set")

    def teardown(self):
        DB.drop_database('atlas-test')


class TestAddNewVariantSet(BaseTest):

    def test_add_new_vcf_variant_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id)
        vcf.add_to_database()
        assert VariantSet.objects().count() == 1
        vs = VariantSet.objects()[0]
        assert vs.name == "test.vcf"


class TestAddNewVariantSetMetaData(BaseTest):

    def test_add_new_vcf_variant_set(self):
        vcf = VCF(
            f="tests/vcf_tests/test.vcf",
            reference_set_id=self.reference_set.id)
        vcf.add_to_database()
        assert VariantSetMetadata.objects().count() >= 1
        assert VariantSetMetadata.objects(key="KMER").count() == 1
