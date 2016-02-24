import vcf
import os.path
from mongoengine import DoesNotExist
from mongoengine import NotUniqueError
from atlas.variants import CallSet
from atlas.variants import Variant
from atlas.variants import VariantSet
from atlas.variants import VariantSetMetadata
from atlas.variants import Call
from atlas.references import ReferenceSet
from atlas.references import Reference
from atlas.utils import make_var_hash

GLOBAL_VARIANT_SET_NAME = "global_atlas"


class VCF(object):

    def __init__(self, f, reference_set_id, method):
        self.f = f
        self.reference_set = ReferenceSet.objects.get(id=reference_set_id)
        self.references = self._create_reference_lookup()
        self.method = method
        self.vcf_variant_set = None
        self.vcf_reader = None
        self.variants = []
        self.calls = []
        self.call_sets = {}

    def _create_reference_lookup(self):
        refs = {}
        for reference in self.reference_set.references:
            refs[reference.name] = reference
        return refs

    def add_to_database(self):
        with open(self.f, 'r') as infile:
            self.vcf_reader = vcf.Reader(infile)
            self._create_new_variant_set()
            self._create_variant_set_meta_data()
            self._create_call_sets()
            self._create_variants_and_calls()

    def _create_variants_and_calls(self):
        for record in self.vcf_reader:
            if not record.FILTER and self._is_record_valid(record):
                v = self._get_or_create_variant(record)
                for call in record.samples:
                    genotype_likelihoods = self._get_genotype_likelihoods(call)
                    c = Call.create(
                        variant=v,
                        call_set=self.call_sets[call.sample],
                        genotype=call['GT'],
                        genotype_likelihoods=genotype_likelihoods)
                    self.calls.append(c)
        Call.objects.insert(self.calls)

    def _get_or_create_variant(self, record):
        try:
            var_hash = make_var_hash(record.REF, record.POS, [
                str(a) for a in record.ALT])
            return Variant.objects.get(var_hash=var_hash)
        except DoesNotExist:
            try:
                reference = self.references[record.CHROM]
            except KeyError as e:
                raise KeyError(
                    "Reference %s cannot be found in reference set %s (%s). Please add it to the database." %
                    (record.CHROM, self.reference_set.id, self.reference_set.name))

            v = Variant.create_and_save(
                variant_sets=self.variant_sets,
                start=record.POS,
                reference_bases=record.REF,
                alternate_bases=[
                    str(a) for a in record.ALT],
                reference=reference,
                names=[record.ID])
            return v

    def _create_new_variant_set(self):
        self.vcf_variant_set = VariantSet.create_and_save(
            name=os.path.basename(
                self.f),
            reference_set=self.reference_set)

    def _create_call_sets(self):
        for sample in self.vcf_reader.samples:
            try:
                cs = CallSet.create_and_save(
                    name="_".join(
                        [
                            sample,
                            self.method]),
                    variant_sets=self.variant_sets,
                    sample_id=sample,
                    info={
                        "variant_caller": self.method})
            except NotUniqueError:
                raise ValueError(
                    "There is already a call set for sample %s with method %s " %
                    (sample, self.method))
            else:
                self.call_sets[sample] = cs

    @property
    def variant_sets(self):
        return [self.vcf_variant_set, self.global_variant_set]

    @property
    def global_variant_set(self):
        try:
            vs = VariantSet.objects.get(name=GLOBAL_VARIANT_SET_NAME)
        except:
            vs = VariantSet.create_and_save(
                name=GLOBAL_VARIANT_SET_NAME,
                reference_set=self.reference_set)
        return vs

    def _create_variant_set_meta_data(self):
        for variant_set in self.variant_sets:
            for k, v in self.vcf_reader.metadata.items():
                if not VariantSetMetadata.objects(
                        key=k,
                        variant_set=variant_set):
                    vsm = VariantSetMetadata.create_and_save(
                        key=k,
                        value="metadata",
                        type="metadata",
                        variant_set=variant_set)
            for k, v in self.vcf_reader.infos.items():
                if not VariantSetMetadata.objects(
                        key=k,
                        variant_set=variant_set):
                    vsm = VariantSetMetadata.create_and_save(
                        key=k,
                        value="infos",
                        type=v.type,
                        variant_set=variant_set,
                        number=int(
                            v.num),
                        description=v.desc)
            for k, v in self.vcf_reader.filters.items():
                if not VariantSetMetadata.objects(
                        key=k,
                        variant_set=variant_set):
                    vsm = VariantSetMetadata.create_and_save(
                        key=k,
                        value="filters",
                        type="filters",
                        variant_set=variant_set,
                        description=v.desc)
            for k, v in self.vcf_reader.formats.items():
                if not VariantSetMetadata.objects(
                        key=k,
                        variant_set=variant_set):
                    vsm = VariantSetMetadata.create_and_save(
                        key=k,
                        value="formats",
                        type=v.type,
                        variant_set=variant_set,
                        number=int(
                            v.num),
                        description=v.desc)

    def _is_record_valid(self, record):
        valid = True
        for sample in record.samples:
            if sample["GT"] is None:
                valid = False
            else:
                if sum([int(i) for i in sample['GT'].split('/')]) < 2:
                    valid = False
            try:
                if sample["GT_CONF"] <= 1:
                    valid = False
            except AttributeError:
                pass
        return valid

    def _get_genotype_likelihoods(self, sample):
        try:
            genotype_likelihoods = [float(i) for i in sample['GL'].split(",")]
        except:
            genotype_likelihoods = [0, 0, 0]
            genotype_likelihoods[
                sum([int(i) for i in sample['GT'].split('/')])] = sample["GT_CONF"]
        return genotype_likelihoods
