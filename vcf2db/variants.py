import datetime
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField
from mongoengine import FloatField


# class VariantSetMetadata(Document)
#  {
#   /** The top-level key. */
#   string key;

#   /** The value field for simple metadata. */
#   string value;

#   /**
#   User-provided ID field, not enforced by this API.
#   Two or more pieces of structured metadata with identical
#   id and key fields are considered equivalent.
#   */
#   string id;

#   /** The type of data. */
#   string type;

#   /**
#   The number of values that can be included in a field described by this
#   metadata.
#   */
#   string number;

#   /** A textual description of this metadata. */
#   string description;

#   /** Remaining structured metadata key-value pairs. */
#   map<array<string>> info = {};
# }

"""
`Variant` and `CallSet` both belong to a `VariantSet`.
`VariantSet` belongs to a `Dataset`.
The variant set is equivalent to a VCF file.
"""
# record VariantSet {
#   /** The variant set ID. */
#   string id;

#   /** The ID of the dataset this variant set belongs to. */
#   string datasetId;

#   /**
#   The reference set the variants in this variant set are using.
#   */
#   string referenceSetId;

#   /**
#   The metadata associated with this variant set. This is equivalent to
#   the VCF header information not already presented in first class fields.
#   */
#   array<VariantSetMetadata> metadata = [];
# }

# /**
# A `CallSet` is a collection of variant calls for a particular sample.
# It belongs to a `VariantSet`. This is equivalent to one column in VCF.
# */
# record CallSet {

#   /** The call set ID. */
#   string id;

#   /** The call set name. */
#   union { null, string } name = null;

#   /** The sample this call set's data was generated from. */
#   union { null, string } sampleId;

#   /** The IDs of the variant sets this call set has calls in. */
#   array<string> variantSetIds = [];

#   /** The date this call set was created in milliseconds from the epoch. */
#   union { null, long } created = null;

#   /**
#   The time at which this call set was last updated in
#   milliseconds from the epoch.
#   */
#   union { null, long } updated = null;

#   /**
#   A map of additional call set information.
#   */
#   map<array<string>> info = {};
# }


class Call(Document):

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
    call_set_name = StringField()
    genotype = ListField(IntField())
    genotype_likelihood = FloatField()
    variant = ReferenceField('Variant')

    @classmethod
    def create(cls, variant, call_set_name, genotype, genotype_likelihood ):
        if type(genotype) is str:
            genotype = [int(g) for g in genotype.split('/')]
        c = cls(variant = variant, call_set_name = call_set_name, genotype = genotype,
                 genotype_likelihood = genotype_likelihood )
        return c.save()    


class Variant(Document):
    """A `Variant` represents a change in DNA sequence relative to some reference. For example, a variant could represent a SNP or an insertion.
      Variants belong to a `VariantSet`.This is equivalent to a row in VCF."""
    var_id = ReferenceField('VariantSet')
    created_at = DateTimeField(default=datetime.datetime.now)
    updated_at = DateTimeField(default=datetime.datetime.now)
    reference = ReferenceField('Reference')
    start = IntField(required = True) #(0-based)
    end = IntField(required = True) #  The end position (exclusive), resulting in [start, end) closed-open interval.
    reference_bases = StringField()
    alternate_bases = ListField(StringField())
   

    @classmethod
    def create(cls, start = 0, end = 1, reference_bases = "A", alternate_bases = ["T"]):
        v = cls(start = start, end = end, reference_bases = reference_bases , alternate_bases = alternate_bases)
        return v.save()

    @property 
    def alt(self):
        if len(self.alternate_bases) == 1:
            return "".join(self.alternate_bases)
        else:
            return self.alternate_bases

    @property 
    def call(self):
        return Call.objects.get(variant = self)



