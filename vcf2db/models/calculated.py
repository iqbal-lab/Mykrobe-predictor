import datetime
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import ListField
from mongoengine import IntField
from mongoengine import FloatField

class UniqueVariants(Document):
    names = ListField(StringField(), required = True)
    created_at = DateTimeField(required = True, default=datetime.datetime.now)

    @classmethod
    def create(cls, names, reference = None):
        c = cls(names = names)
        return c.save()

class VariantFreq(Document):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['name']
                }                                                   
                ]
            }    
    name = StringField()
    count = IntField()
    total_samples = IntField()
    freq = FloatField()
    created_at = DateTimeField(required = True, default=datetime.datetime.now)
    
    start = IntField()
    reference_bases = StringField(required = False)
    alternate_bases = ListField(StringField(), required = False)    


    @classmethod
    def create(cls, name, count, total_samples, reference_bases, start, alternate_bases):
        freq = float(count) / float(total_samples)
        return cls(name = name,
                   count = count,
                   total_samples = total_samples,
                   freq = freq, 
                   start = int(start),
                   reference_bases = reference_bases, 
                   alternate_bases = alternate_bases
                   ).save()

    def __str__(self):
        return "".join([self.reference_bases, str(self.start), "/".join(self.alternate_bases)])

    def __repr__(self):
        return "".join([self.reference_bases, str(self.start), "/".join(self.alternate_bases)])