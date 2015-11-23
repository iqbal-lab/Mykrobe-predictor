import datetime
from atlas.vcf2db.models.base import make_hash
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import ListField
from mongoengine import IntField
from mongoengine import FloatField
from mongoengine import ReferenceField

class UniqueVariants(Document):
    names = ListField(StringField(), required = True)
    created_at = DateTimeField(required = True, default=datetime.datetime.now)

    @classmethod
    def create(cls, names, reference = None):
        c = cls(names = names)
        return c.save()
        
class VariantPanel(Document):

    meta = {'indexes': [
                {
                    'fields' : ['variant']
                },
                {
                    'fields' : ['name_hash']
                }                                                                  
                ]
            } 

    variant = ReferenceField('VariantFreq')
    ref = StringField()    
    alts = ListField(StringField())
    name = StringField(default = None)
    name_hash = StringField(unique = True)

    @classmethod
    def create(cls, variant, vf, ref, alts):
        return cls(
        variant = vf, 
        ref = ref, 
        alts = alts,
        name = variant.name,
        name_hash = variant.name_hash
        )

    def update(self, ref, alts):
        self.ref = ref
        self.alts = alts
        return self.save()
        
    def __repr__(self):
        return "PANAL %s" % self._name

class VariantFreq(Document):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['name_hash']
                }                                                   
                ]
            }    
    name = StringField()
    name_hash = StringField()
    # count = IntField()
    # total_samples = IntField()
    freq = FloatField()
    created_at = DateTimeField(required = True, default=datetime.datetime.now)
    
    start = IntField()
    reference_bases = StringField(required = False)
    alternate_bases = ListField(StringField(), required = False)    


    @classmethod
    def create(cls, name, reference_bases, start, alternate_bases, name_hash = None, total_samples= None):
        if name_hash is None:
            name_hash = make_hash(name)
        return cls(name = name,
                    name_hash = name_hash,
                   # count = count,
                   # total_samples = total_samples,
                   # freq = freq, 
                   start = int(start),
                   reference_bases = reference_bases, 
                   alternate_bases = alternate_bases
                   )

    def __str__(self):
        return "".join([self.reference_bases, str(self.start), "/".join(self.alternate_bases)])

    def __repr__(self):
        return "".join([self.reference_bases, str(self.start), "/".join(self.alternate_bases)])