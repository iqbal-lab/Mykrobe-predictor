#! /usr/bin/env python
from mongoengine import connect
connect('atlas')
from models import Variant
from models import UniqueVariants

data = Variant.objects().distinct('name')
UniqueVariants.create(names = data)