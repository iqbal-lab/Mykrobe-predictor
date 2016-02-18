import datetime
import json


from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField
from mongoengine import FloatField
from mongoengine import DictField
from mongoengine import GenericReferenceField

from atlas.variants.models.base import CreateAndSaveMixin
from atlas.utils import split_var_name
from atlas.utils import make_var_hash

# Based on ga4gh Variant schema http://ga4gh.org/#/schemas feb 2016 with
# ocassional changes


class VariantSetMetadata(Document, CreateAndSaveMixin):
    key = StringField()
    value = StringField()
    type = StringField()
    # The number of values that can be included in a field described by this
    # metadata.
    number = IntField()
    description = StringField()
    info = DictField()
    variant_set = ReferenceField("VariantSet")

    @classmethod
    def create(
            cls,
            key,
            value,
            type,
            variant_set,
            number=None,
            description=None,
            info=None):
        return cls(key=key, value=value, type=type, description=description,
                   info=info, variant_set=variant_set)


class VariantSet(Document, CreateAndSaveMixin):

    """
    `Variant` and `CallSet` both belong to a `VariantSet`.
    `VariantSet` belongs to a `Dataset`.
    The variant set is equivalent to a VCF file.
    """
    name = StringField(required=True, unique=True)
    dataset = ReferenceField('Dataset')
    reference_set = ReferenceField('ReferenceSet', required=True)

    @classmethod
    def create(cls, name, reference_set=None, dataset=None):
        c = cls(name=name, reference_set=reference_set, dataset=dataset)
        return c.save()

    # Optional metadata associated with this variant set.
    # This array can be used to store information about the variant set, such as information found
    # in VCF header fields, that isn't already available in first class fields
    # such as "name".
    @property
    def metadata(self):
        return VariantSetMetadata.objects(variant_set=self)


class CallSet(Document, CreateAndSaveMixin):

    """
         A `CallSet` is a collection of variant calls for a particular sample.
         It belongs to a `VariantSet`. This is simillar to one column in VCF.


    """
    name = StringField(required=True, default=None, unique = True)
    sample_id = StringField(required=True)
    created_at = DateTimeField(default=datetime.datetime.now)
    updated_at = DateTimeField(required=True, default=datetime.datetime.now)
    # Break from ga4gh schema -
    variant_sets = ListField(ReferenceField('VariantSet'))
    # When can a call set exist in multiple variant sets? If you have a set of
    # calls that you want to add to multiple variant sets. I think this demands
    # that a variant can exist in multiple variant sets, something not allowed
    # by ga4gh schema.
    info = DictField()

    @classmethod
    def create(cls, name, variant_sets, sample_id=None, info={}):
        c = cls(name=name, variant_sets=variant_sets, sample_id=sample_id,
                info=info)
        return c.save()


def convert_string_gt_to_list_int_gt(variant, genotype):
    return [int(i) for i in genotype.split('/')]


class Call(Document, CreateAndSaveMixin):
    meta = {'indexes': [
        {
            'fields': ['call_set']
        },
        {
            'fields': ['variant']
        }
    ]
    }
    """
    A `Call` represents the determination of genotype with respect to a
    particular `Variant`.
    It may include associated information such as quality
    and phasing. For example, a call might assign a probability of 0.32 to
    the occurrence of a SNP named rs1234 in a call set with the name NA12345.
    """
    """
    The name of the call set this variant call belongs to.
    If this field is not present, the ordering of the call sets from a
    `SearchCallSetsRequest` over this `VariantSet` is guaranteed to match
    the ordering of the calls on this `Variant`.
    The number of results will also be the same.
    """
    variant = ReferenceField('Variant', required=True)  # Not in ga4gh
    call_set = ReferenceField('CallSet', required=True)
    """
    The genotype of this variant call.

    A 0 value represents the reference allele of the associated `Variant`. Any
    other value is a 1-based index into the alternate alleles of the associated
    `Variant`.

    If a variant had a referenceBases field of "T", an alternateBases
    value of ["A", "C"], and the genotype was [2, 1], that would mean the call
    represented the heterozygous value "CA" for this variant. If the genotype
    was instead [0, 1] the represented value would be "TA". Ordering of the
    genotype values is important if the phaseset field is present.
    """

    genotype = ListField(IntField())
    """
    The genotype likelihoods for this variant call. Each array entry
    represents how likely a specific genotype is for this call as
    log10(P(data | genotype)), analogous to the GL tag in the VCF spec. The
    value ordering is defined by the GL tag in the VCF spec (below)

    GL : genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible
    genotypes given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same
    ploidy is expected and the canonical order is used; without GT field, diploidy is assumed. If A is the allele in
    REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by:
    F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the
    ordering is: AA,AB,BB,AC,BC,CC, etc. For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Floats)

    """
    genotype_likelihoods = ListField(FloatField())
    # If this field is not null, this variant call's genotype ordering implies
    # the phase of the bases and is consistent with any other variant calls on
    # the same contig which have the same phaseset string.
    phaseset = GenericReferenceField(default=None)
    info = DictField()

    @classmethod
    def create(cls, variant, call_set, genotype, genotype_likelihoods=[],
               phaseset=None, info={}):
        if isinstance(genotype, str):
            genotype = convert_string_gt_to_list_int_gt(variant, genotype)
        cls._check_genotype_likelihood_length(genotype_likelihoods, variant)
        return cls(
            variant=variant,
            call_set=call_set,
            genotype=genotype,
            genotype_likelihoods=genotype_likelihoods,
            phaseset=phaseset,
            info=info)

    @classmethod
    def _check_genotype_likelihood_length(cls, genotype_likelihood, variant):
        if len(variant.alternate_bases) == 1:
            if not len(genotype_likelihood) == 3:
                raise ValueError(
                    "Biallelic sites should have 3 genotype likelihoods. AA,AB,BB")
        elif len(variant.alternate_bases) == 2:
            if not len(genotype_likelihood) == 6:
                raise ValueError(
                    "Biallelic sites should have 6 genotype likelihoods. AA,AB,BB,AC,BC,CC, etc")
        else:
            raise NotImplementedError(
                "Haven't implemented check for > triallelic sites")

    @property
    def call_set_name(self):
        return self.call_set.name


def lazyprop(fn):
    attr_name = '_' + fn.__name__

    @property
    def _lazyprop(self):
        if not getattr(self, attr_name):
            setattr(self, attr_name, fn(self))
            self.save()
        return getattr(self, attr_name)
    return _lazyprop


class Variant(Document, CreateAndSaveMixin):
    meta = {'indexes': [
        {
            'fields': ['start']
        },
        {
            'fields': ['var_hash']
        },
        {
            'fields': ['variant_sets']
        }
    ]
    }
    """A `Variant` represents a change in DNA sequence relative to some reference.
       For example, a variant could represent a SNP or an insertion.
      Variants belong to a `VariantSet`. This is simillar to a row in VCF.

      However, breaking from ga4gh we're allowing a variant belong to multiple
      VariantSets to allow a Variant to belong to multiple "Sets" of variants.
      """
    # Here, we've broken from ga4gh as they demand every variant belongs to a
    # single variant set. See CallSet.
    variant_sets = ListField(ReferenceField('VariantSet'), required=True)
    # The var_hash is a unique description on a variant. We use the hash of
    # "ref+pos+alt".
    var_hash = StringField(required=True)
    names = ListField(StringField())
    created_at = DateTimeField(required=True, default=datetime.datetime.now)
    updated_at = DateTimeField(required=True, default=datetime.datetime.now)
    start = IntField(required=True)  # (0-based)
    # The end position (exclusive), resulting in [start, end) closed-open
    # interval.
    end = IntField(required=False)
    reference_bases = StringField(required=True)
    alternate_bases = ListField(StringField(), required=True)
    info = DictField()

    # Each variant is defined against a single reference.
    # We can't have the same variant more than once i
    reference = ReferenceField(
        "Reference",
        required=True,
        unique_with="var_hash")

    _length = IntField(required=False, default=None)

    @classmethod
    def create(cls, variant_sets, start, reference_bases,
               alternate_bases, reference, end=None,
               names=[]):
        
        return cls(variant_sets=variant_sets,
                   start=start,
                   end=end,
                   reference_bases=reference_bases,
                   alternate_bases=alternate_bases,
                   reference=reference,
                   var_hash=make_var_hash(reference_bases, start, alternate_bases))

    @property
    def calls(self):
        # The variant calls for this particular variant. Each one represents the
        # determination of genotype with respect to this variant. `Call`s in this array
        # are implicitly associated with this `Variant`.
        return Call.objects(variant=self)

    @property
    def var_name(self):
        return "".join(
            [reference_bases, str(start), "/".join(alternate_bases)])

    @lazyprop
    def length(self):
        return abs(len(self.reference_bases) -
                   max([len(a) for a in self.alternate_bases]))

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return self.name_hash == other.name_hash

    @property
    def alt(self):
        if len(self.alternate_bases) == 1:
            return "".join(self.alternate_bases)
        else:
            return self.alternate_bases

    @property
    def is_indel(self):
        """ Return whether or not the variant is an INDEL """
        if len(self.reference_bases) > 1:
            return True
        for alt in self.alternate_bases:
            if alt is None:
                return True
            elif len(alt) != len(self.reference_bases):
                return True
        return False

    @property
    def is_deletion(self):
        """ Return whether or not the INDEL is a deletion """
        # if multiple alts, it is unclear if we have a transition
        if len(self.alternate_bases) > 1:
            return False

        if self.is_indel:
            # just one alt allele
            alt_allele = self.alternate_bases[0]
            if alt_allele is None:
                return True
            if len(self.reference_bases) > len(alt_allele):
                return True
            else:
                return False
        else:
            return False

    @property
    def is_insertion(self):
        if len(self.alternate_bases) > 1:
            return False
        if self.is_indel:
            # just one alt allele
            alt_allele = self.alternate_bases[0]
            if alt_allele is None:
                return False
            if len(alt_allele) > len(self.reference_bases):
                return True
            else:
                return False
        else:
            return False
