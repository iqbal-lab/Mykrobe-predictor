import vcf
import os.path
from atlas.variants import CallSet
from atlas.variants import Reference
from atlas.variants import Variant
from atlas.variants import VariantSet
from atlas.variants import VariantSetMetadata
from atlas.variants import Call
from atlas.references import ReferenceSet


class VCF(object):

    def __init__(self, f, reference_set_id):
        self.f = f
        self.reference_set = ReferenceSet.objects.get(id=reference_set_id)
        self.variant_set = None
        self.vcf_reader = None

    def add_to_database(self):
        with open(self.f, 'r') as infile:
            self.vcf_reader = vcf.Reader(infile)
            self._create_new_variant_set()
            self._create_variant_set_meta_data()

    def _create_new_variant_set(self):
        self.variant_set = VariantSet.create_and_save(
            name=os.path.basename(
                self.f),
            reference_set=self.reference_set)

    def _create_variant_set_meta_data(self):
        variant_set_metadatas = []
        for k, v in self.vcf_reader.metadata.items():
            vsm = VariantSetMetadata.create(key=k,
                                            value="metadata",
                                            type="metadata",
                                            variant_set=self.variant_set)
            variant_set_metadatas.append(vsm)
        for k, v in self.vcf_reader.infos.items():
            vsm = VariantSetMetadata.create(key=k,
                                            value="infos",
                                            type=v.type,
                                            variant_set=self.variant_set,
                                            number=int(v.num),
                                            description=v.desc)
            variant_set_metadatas.append(vsm)
        for k, v in self.vcf_reader.filters.items():
            vsm = VariantSetMetadata.create(key=k,
                                            value="filters",
                                            type="filters",
                                            variant_set=self.variant_set,
                                            description=v.desc)
            variant_set_metadatas.append(vsm)
        for k, v in self.vcf_reader.formats.items():
            vsm = VariantSetMetadata.create(key=k,
                                            value="formats",
                                            type=v.type,
                                            variant_set=self.variant_set,
                                            number=int(v.num),
                                            description=v.desc)
            variant_set_metadatas.append(vsm)
        VariantSetMetadata.objects.insert(variant_set_metadatas)
