from mongoengine import Document
from mongoengine import StringField
from mongoengine import IntField
from mongoengine import FloatField
from mongoengine import ListField
from mongoengine import ReferenceField

class TypedPresence(Document):


  name = StringField()
  version = IntField()
  percent_coverage = FloatField()
  median_depth = StringField()
  alt_names = ListField(StringField())
  name_hash = StringField(unique_with = "call_set")
  call_set = ReferenceField('CallSet')

  @classmethod
  def create_object(cls, name, version, percent_coverage, median_depth, alt_names = []):
    if not alt_names:
      alt_names = [self.version]
    return cls(name = name,
      version = version,
      percent_coverage = percent_coverage,
      median_depth = median_depth,
      alt_names = alt_names)      

  def __str__(self):
      return self.version

  def __repr__(self):
      return self.version 

  @property 
  def gene_version(self):
      return "-".join([self.name, self.version])

  def to_dict(self):
      d  = {  "name" : self.name,
              "alt_name" : ",".join(self.alt_names),
              "covg" : {"percent_coverage" : self.percent_coverage, 
                        "median_depth" : self.median_depth
                        }
            }
      return d        