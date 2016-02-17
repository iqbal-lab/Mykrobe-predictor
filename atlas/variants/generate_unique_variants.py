#! /usr/bin/env python
from models import Variant
from models import VariantFreq
from models import VariantSet
from mongoengine import connect

connect('atlas')
total_samples = VariantSet.objects.count()
results = Variant.objects().aggregate(
    {
        "$group": {
            "_id": {
                "name": "$name"}, "count": {
                    "$sum": 1}}}, {
                        "$match": {
                            "count": {
                                "$gt": 1}}})
for result in results:
    var = Variant.objects(name=result["_id"]["name"])[0]
    VariantFreq.create(name=result["_id"]["name"],
                       count=result["count"],
                       total_samples=total_samples,
                       start=var['start'],
                       reference_bases=var['reference_bases'],
                       alternate_bases=var['alternate_bases']
                       )
