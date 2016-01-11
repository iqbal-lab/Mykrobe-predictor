from mongoengine import Document
from mongoengine import StringField
from mongoengine import IntField
from mongoengine import FloatField
from mongoengine import ListField
from mongoengine import ReferenceField

class SequenceCoverage(Document):

  "Summary of kmer coverage of sequence. e.g output of color covearges"

  name = StringField()
  version = IntField()
  percent_coverage = FloatField()
  median_depth = StringField()
  alt_names = ListField(StringField())
  percent_coverage_threshold = IntField(default = 30)
  gt = StringField()
  induced_resistance = ListField(StringField())
  length = IntField()
  copy_number = FloatField()


  @classmethod
  def create_object(cls, name, percent_coverage,
                   median_depth, version = 1, alt_names = [],
                   percent_coverage_threshold = 30,
                   length = None):
    if not alt_names:
        alt_names = ["-".join([name, str(version)])]
    return cls(name = name,
      version = version,
      percent_coverage = percent_coverage,
      median_depth = median_depth,
      alt_names = alt_names,
      percent_coverage_threshold = percent_coverage_threshold,
      length = length)      

  def __str__(self):
      return str(self.to_dict())

  # def __repr__(self):
  #     return str(self.to_dict())

  @property 
  def gene_version(self):
      return "-".join([self.name, self.version])

  def to_dict(self):
      d  = {  "name" : self.name,
              "alt_name" : ",".join(self.alt_names),
              "gt" : self.gt,
              "covg" : {"percent_coverage" : self.percent_coverage, 
                        "median_depth" : self.median_depth
                        },
              "induced_resistance" : self.induced_resistance
            }
      return d 

  def add_induced_resistance(self, drug):
      if drug not in self.induced_resistance:
          self.induced_resistance.append(drug)

  def set_genotype(self, gt):
      self.gt = gt

  def set_copy_number(self, cn):
      self.copy_number = cn

  @property
  def alt_name(self):
      return "".join(self.alt_names)
