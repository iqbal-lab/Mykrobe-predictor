import datetime
import json

from mongoengine import Document
from mongoengine import StringField
from mongoengine import DateTimeField
from mongoengine import IntField
from mongoengine import ReferenceField
from mongoengine import ListField
from mongoengine import FloatField

from atlas.utils import split_var_name
from atlas.utils import make_hash


# class TypedVariant(Document):
#     meta = {'indexes': [
#         {
#             'fields': ['start']
#         },
#         {
#             'fields': ['name_hash']
#         },
#         {
#             'fields': ['call_set']
#         }
#     ]
#     }
#     name = StringField()
#     alt_name = StringField()
#     name_hash = StringField(unique_with="call_set")
#     reference_median_depth = IntField()
#     alternate_median_depth = IntField()
#     reference_percent_coverage = FloatField()
#     alternate_percent_coverage = FloatField()
#     created_at = DateTimeField(required=True, default=datetime.datetime.now)

#     start = IntField()
#     alt_index = IntField()
#     reference_bases = StringField()
#     alternate_bases = StringField()
#     call_set = ReferenceField('CallSet')
#     gt = StringField()
#     copy_number = FloatField()
#     induced_resistance = ListField(StringField)

#     @classmethod
#     def create_object(cls,
#                       reference_percent_coverage,
#                       alternate_percent_coverage,
#                       reference_median_depth,
#                       alternate_median_depth,
#                       name="", call_set=None, gt=None,
#                       alt_name=None, alt_index=None):
#         reference_bases, start, alternate_bases = split_var_name(name)
#         if reference_median_depth is None:
#             reference_median_depth = 0
#         if alternate_median_depth is None:
#             alternate_median_depth = 0

#         return cls(name=name,
#                    name_hash=make_hash(name),
#                    reference_bases=reference_bases,
#                    start=start,
#                    alternate_bases=alternate_bases,
#                    call_set=call_set,
#                    reference_percent_coverage=int(reference_percent_coverage),
#                    alternate_percent_coverage=int(alternate_percent_coverage),
#                    reference_median_depth=int(reference_median_depth),
#                    alternate_median_depth=int(alternate_median_depth),
#                    gt=gt,
#                    alt_name=alt_name,
#                    alt_index=alt_index
#                    )

#     @classmethod
#     def create(cls, **kwargs):
#         return cls.create_object(**kwargs).save()

#     def set_genotype(self, gt):
#         self.gt = gt

#     def set_copy_number(self, cn):
#         self.copy_number = cn

#     def add_induced_resistance(self, drug):
#         if drug not in self.induced_resistance:
#             self.induced_resistance.append(drug)

#     def to_dict(self):
#         d = {
#             "name": self.name,
#             "alt_name": self.alt_name,
#             "gt": self.gt,
#             "covg": {
#                 "reference_percent_coverage": self.reference_percent_coverage,
#                 "alternate_percent_coverage": self.alternate_percent_coverage,
#                 "reference_median_depth": self.reference_median_depth,
#                 "alternate_median_depth": self.alternate_median_depth},
#             "induced_resistance": self.induced_resistance}
#         return d
