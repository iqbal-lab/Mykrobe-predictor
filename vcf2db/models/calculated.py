import datetime
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import ListField

class UniqueVariants(Document):
    """A `Variant` represents a change in DNA sequence relative to some reference. For example, a variant could represent a SNP or an insertion.
      Variants belong to a `VariantSet`.This is equivalent to a row in VCF."""
    names = ListField(StringField(), required = True)
    created_at = DateTimeField(required = True, default=datetime.datetime.now)

    @classmethod
    def create(cls, names, reference = None):
        c = cls(names = names)
        return c.save()