import datetime
from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField
from mongoengine import FloatField

from atlas.utils import split_var_name
from atlas.utils import make_hash

class GenotypedVariant(Document):
    meta = {'indexes': [
                {
                    'fields' : ['start']
                },
                {
                    'fields' : ['name_hash']
                },
                {
                    'fields' : ['call_set']
                }                                                  
                ]
            }    
    name = StringField()
    alt_name = StringField()
    name_hash = StringField(unique_with = "call_set")
    ref_coverage = IntField()
    alt_coverage = IntField()
    ref_pnz = FloatField()
    alt_pnz = FloatField()    
    created_at = DateTimeField(required = True, default = datetime.datetime.now)
    
    start = IntField()
    reference_bases = StringField()
    alternate_bases = StringField()
    call_set = ReferenceField('CallSet')
    gt = StringField()

    @classmethod
    def create_object(cls, name, call_set, ref_pnz, alt_pnz, ref_coverage, alt_coverage, gt, alt_name = None):
        reference_bases, start, alternate_bases = split_var_name(name)
        if ref_coverage is None:
            ref_coverage = 0
        if alt_coverage is None:
            alt_coverage = 0   

        return cls( name = name, 
                    name_hash = make_hash(name),
                    reference_bases = reference_bases, 
                    start = start, 
                    alternate_bases = alternate_bases,
                    call_set = call_set,
                    ref_pnz = int(ref_pnz), 
                    alt_pnz = int(alt_pnz),
                    ref_coverage = int(ref_coverage),
                    alt_coverage = int(alt_coverage),
                    gt = gt,
                    alt_name = alt_name           
                    )   

    @classmethod
    def create(cls, name, call_set, ref_pnz, alt_pnz, ref_coverage, alt_coverage, gt):
        return cls.create_object(name, call_set, ref_pnz, alt_pnz, ref_coverage, alt_coverage, gt).save()

    def to_dict(self):
        d  = {  "name" : self.name,
                "alt_name" : self.alt_name,
                "gt" : self.gt,
                "covg" : {"reference_percent_coverage" : self.ref_pnz, 
                          "alternate_percent_coverage" : self.alt_pnz, 
                          "reference_median_depth": self.ref_coverage,
                          "alternate_median_depth" : self.alt_coverage
                          }
              }
        return d