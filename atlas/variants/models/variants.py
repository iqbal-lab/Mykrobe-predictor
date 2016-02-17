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
from atlas.utils import make_hash
from atlas.utils import split_var_name

## Based on ga4gh Variant schema http://ga4gh.org/#/schemas feb 2016 with ocassional changes

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
    def create(cls, name, dataset_id = None, reference_set_id= None):
        c = cls(name = name)
        return c.save()    

class VariantSet(Document, CreateAndSaveMixin):
    """
    `Variant` and `CallSet` both belong to a `VariantSet`.
    `VariantSet` belongs to a `Dataset`.
    The variant set is equivalent to a VCF file.
    """
    name = StringField(required = True, unique = True)
    dataset = ReferenceField('Dataset')
    reference = ReferenceField('ReferenceSet')

    @classmethod
    def create(cls, name, dataset = None, reference_set= None):
        c = cls(name = name)
        return c.save() 

    # Optional metadata associated with this variant set.
    # This array can be used to store information about the variant set, such as information found
    # in VCF header fields, that isn't already available in first class fields such as "name".
    @property
    def metadata(self):
        return VariantSetMetadata.objects(variant_set = self)


class CallSet(Document, CreateAndSaveMixin):
    """
         A `CallSet` is a collection of variant calls for a particular sample.
         It belongs to a `VariantSet`. This is simillar to one column in VCF.


    """
    name = StringField(required = True, default = None)
    sample_id = StringField(required = True)
    created_at = DateTimeField(default = datetime.datetime.now)
    updated_at = DateTimeField(required = True, default=datetime.datetime.now)
    variant_sets = ListField(ReferenceField('VariantSet')) ## Break from ga4gh schema -
    ## When can a call set exist in multiple variant sets? If you have a set of
    ## calls that you want to add to multiple variant sets. I think this demands
    ## that a variant can exist in multiple variant sets, something not allowed by ga4gh schema.

    info = DictField()

    @classmethod
    def create(cls, name, sample_id = None):
        c = cls(name = name, sample_id = sample_id)
        return c.save() 

class Call(Document, CreateAndSaveMixin):
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
    genotype_likelihood = ListField(FloatField())
    # If this field is not null, this variant call's genotype ordering implies
    # the phase of the bases and is consistent with any other variant calls on
    # the same contig which have the same phaseset string.
    phaseset = GenericReferenceField(default = None)
    info = DictField()
    variant = ReferenceField('Variant', required = True) # Not in ga4gh


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

class Variant(Document, CreateAndSaveMixin):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['var_hash']
                },
                {
                    'fields' : ['variant_sets']
                }                                                               
                ]
            }    
    """A `Variant` represents a change in DNA sequence relative to some reference. 
       For example, a variant could represent a SNP or an insertion.
      Variants belong to a `VariantSet`. This is simillar to a row in VCF.

      However, breaking from ga4gh we're allowing a variant belong to multiple 
      VariantSets to allow a Variant to belong to multiple "Sets" of variants.
      """
    ## Here, we've broken from ga4gh as they demand every variant belongs to a 
    ## single variant set. See CallSet.       
    variant_sets = ListField(ReferenceField('VariantSet'), required = True) 
    ## The var_hash is a unique description on a variant. We use the hash of "ref+pos+alt".
    var_hash = StringField(required = True)
    names = ListField(StringField())
    created_at = DateTimeField(required = True, default=datetime.datetime.now)
    updated_at = DateTimeField(required = True, default=datetime.datetime.now)
    start = IntField(required = True) #(0-based)
    end = IntField(required = False) #  The end position (exclusive), resulting in [start, end) closed-open interval.
    reference_bases = StringField(required = True)
    alternate_bases = ListField(StringField(), required = True)
    info = DictField()

    ## Each variant is defined against a single reference. 
    ## We can't have the same variant more than once i
    reference = ReferenceField("Reference", required = True, unique_with = "var_hash")

    _length = IntField(required = False, default = None)  

    @classmethod
    def create(cls, variant_sets, start,  reference_bases,
                      alternate_bases, reference, end = None, names = []):
        name = "".join([reference_bases,str(start),"/".join(alternate_bases)])
        return cls(variant_sets = variant_sets,
                   start = start,
                   end = end,
                   reference_bases = reference_bases,
                   alternate_bases = alternate_bases,
                   reference = reference,
                   var_hash = make_hash(name))

    @property
    def calls(self):
      # The variant calls for this particular variant. Each one represents the
      # determination of genotype with respect to this variant. `Call`s in this array
      # are implicitly associated with this `Variant`.
       return Call.objects(variant = self)        

    # @classmethod
    # def create(cls, variant_set, start,  reference_bases, alternate_bases,
    #                  reference, end = None):
    #     return cls().create_object(variant_set = variant_set,
    #                         start = start,
    #                         reference_bases = reference_bases,
    #                         alternate_bases = alternate_bases,
    #                         reference = reference,
    #                         end = end).save()

    @property
    def long_name(self):
        return "".join([reference_bases,str(start),"/".join(alternate_bases)])


    @lazyprop 
    def length(self):
        return abs(len(self.reference_bases) - max([len(a) for a in self.alternate_bases]))

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
