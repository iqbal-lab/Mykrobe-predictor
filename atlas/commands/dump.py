#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import datetime
import logging
LOGGER = logging.getLogger("logger")


from mongoengine import connect

from atlas.schema import Variant
from atlas.panelgeneration import AlleleGenerator
from atlas.panelgeneration import make_variant_probe


def run(parser, args):
    connect('atlas-%s' % (args.db_name))
    if args.verbose:
        LOGGER.setLevel(level=logging.DEBUG)
    else:
        LOGGER.setLevel(level=logging.ERROR)
    al = AlleleGenerator(
        reference_filepath=args.reference_filepath,
        kmer=args.kmer)
    for variant in Variant.snps():
        variant_panel = make_variant_probe(al, variant, args.kmer)
        sys.stdout.write(
            ">ref-%s?num_alts=%i&ref=%s\n" %
            (variant_panel.variant.var_hash, len(
                variant_panel.alts), variant_panel.variant.reference.id))
        sys.stdout.write("%s\n" % variant_panel.ref)
        for a in variant_panel.alts:
            sys.stdout.write(">alt-%s\n" % variant_panel.variant.var_hash)
            sys.stdout.write("%s\n" % a)
