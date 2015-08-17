import datetime
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField

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



