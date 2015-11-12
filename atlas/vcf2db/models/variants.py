import datetime
import re
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField
from mongoengine import FloatField

class GenotypedVariant(Document):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['name']
                },
                {
                    'fields' : ['call_set']
                }                                                  
                ]
            }    
    name = StringField()
    ref_coverage = IntField()
    alt_coverage = IntField()
    ref_pnz = IntField()
    alt_pnz = IntField()    
    created_at = DateTimeField(required = True, default = datetime.datetime.now)
    
    start = IntField()
    reference_bases = StringField()
    alternate_bases = StringField()
    call_set = ReferenceField('CallSet')
    gt = StringField()

    @classmethod
    def create_object(cls, name, call_set, ref_pnz, alt_pnz, ref_coverage, alt_coverage, gt):
        reference_bases, start, alternate_bases = split_var_name(name)
        if ref_coverage is None:
            ref_coverage = 0
        if alt_coverage is None:
            alt_coverage = 0            
        return cls( name = name, 
                    reference_bases = reference_bases, 
                    start = start, 
                    alternate_bases = alternate_bases,
                    call_set = call_set,
                    ref_pnz = int(ref_pnz), 
                    alt_pnz = int(alt_pnz),
                    ref_coverage = int(ref_coverage),
                    alt_coverage = int(alt_coverage),
                    gt = gt                
                    )   

    @classmethod
    def create(cls, name, call_set, ref_pnz, alt_pnz, ref_coverage, alt_coverage, gt):
        return cls.create_object(name, call_set, ref_pnz, alt_pnz, ref_coverage, alt_coverage, gt).save()


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

class VariantSet(Document):
    """
    `Variant` and `CallSet` both belong to a `VariantSet`.
    `VariantSet` belongs to a `Dataset`.
    The variant set is equivalent to a VCF file.
    """
    name = StringField(required = True, unique = True)
    dataset_id = ReferenceField('Dataset')
    reference_set_id = ReferenceField('ReferenceSetId')

    @classmethod
    def create(cls, name, dataset_id = None, reference_set_id= None):
        c = cls(name = name)
        return c.save()   


class CallSet(Document):
    """
         A `CallSet` is a collection of variant calls for a particular sample.
         It belongs to a `VariantSet`. This is equivalent to one column in VCF.
    """
    name = StringField(required = True, unique = True)
    sample_id = StringField()
    created_at = DateTimeField(default = datetime.datetime.now)
    updated_at = DateTimeField(required = True, default=datetime.datetime.now)

#   /** The IDs of the variant sets this call set has calls in. */
#   array<string> variantSetIds = [];

#   map<array<string>> info = {};
# }

    @classmethod
    def create(cls, name, sample_id = None):
        c = cls(name = name, sample_id = sample_id)
        return c.save() 

class Call(Document):
    meta = {'indexes': [
                {
                    'fields' : ['call_set']
                },
                {
                    'fields' : ['variant']
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
    call_set = ReferenceField('CallSet', required = True)
    genotype = ListField(IntField())
    genotype_likelihood = FloatField()
    variant = ReferenceField('Variant', required = True)

    @classmethod
    def create_object(cls, variant, call_set, genotype, genotype_likelihood ):
        if type(genotype) is str:
            genotype = [int(g) for g in genotype.split('/')]
        return {"call_set" : call_set.id,
        "genotype" : genotype,
        "genotype_likelihood" : genotype_likelihood}

    @classmethod
    def create(cls, variant, call_set, genotype, genotype_likelihood ):
        if type(genotype) is str:
            genotype = [int(g) for g in genotype.split('/')]        
        c = cls(variant = variant, call_set = call_set, genotype = genotype,
                 genotype_likelihood = genotype_likelihood ).save()
        return c

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

class Variant(Document):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['name']
                },
                {
                    'fields' : ['variant_set']
                }                                                               
                ]
            }    
    """A `Variant` represents a change in DNA sequence relative to some reference. For example, a variant could represent a SNP or an insertion.
      Variants belong to a `VariantSet`.This is equivalent to a row in VCF."""
    variant_set = ReferenceField('VariantSet', required = True)
    name = StringField(unique_with = "variant_set")
    # created_at = DateTimeField(required = True, default=datetime.datetime.now)
    # updated_at = DateTimeField(required = True, default=datetime.datetime.now)
    reference = ReferenceField('Reference', required = True)
    start = IntField(required = True) #(0-based)
    _length = IntField(required = False, default = None) 
    end = IntField(required = False) #  The end position (exclusive), resulting in [start, end) closed-open interval.
    reference_bases = StringField(required = True)
    alternate_bases = ListField(StringField(), required = True)

    @classmethod
    def create_object(cls, variant_set, start,  reference_bases, alternate_bases, reference, end = None):
        if not len(reference_bases) < 150 or not all([len(a) < 150 for a in alternate_bases ]):
            raise ValueError("INDEL is too large. Atlas can only consider small INDELs < 150 bases for now.")

        name = "".join([reference_bases,str(start),"/".join(alternate_bases)])
        return cls(variant_set = variant_set,
                   start = start, end = end,
                   reference_bases = reference_bases,
                   alternate_bases = alternate_bases,
                   reference = reference,
                   name = name)

    @classmethod
    def create(cls, variant_set, start,  reference_bases, alternate_bases, reference, end = None):
        name = "".join([reference_bases,str(start),"/".join(alternate_bases)])
        return cls(variant_set = variant_set, start = start, end = end, reference_bases = reference_bases,
            alternate_bases = alternate_bases, reference = reference, name = name).save()

    @lazyprop 
    def length(self):
        return abs(len(self.reference_bases) - max([len(a) for a in self.alternate_bases]))

    @property 
    def alt(self):
        if len(self.alternate_bases) == 1:
            return "".join(self.alternate_bases)
        else:
            return self.alternate_bases

    @property 
    def call(self):
        return Call.objects.get(variant = self)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

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


def split_var_name(name):
    items = re.match(r"([A-Z]+)([0-9]+)([A-Z]+)", name, re.I).groups()
    return items[0],int(items[1]),items[2]


