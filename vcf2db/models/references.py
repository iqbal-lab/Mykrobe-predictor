from mongoengine import Document
from mongoengine import StringField
from mongoengine import IntField
from mongoengine import ListField

class Reference(Document):

    length = IntField()
    name = StringField()
    source_accessions = ListField(StringField())

    @classmethod
    def create(cls, name, length, source_accessions ):
    	if type(source_accessions) is str:
    		source_accessions = [source_accessions]
    	return cls(name = name, length = length, source_accessions = source_accessions).save()