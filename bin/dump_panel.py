#! /usr/bin/env python
import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import datetime
import logging
LOGGER = logging.getLogger("logger")
from collections import Counter

from Bio import SeqIO

from mongoengine import connect
import multiprocessing


from atlas.schema import Variant
from atlas.schema import VariantSet
from atlas.utils import split_var_name
from atlas.utils import flatten
from atlas.utils import unique
from atlas.panelgeneration import AlleleGenerator
# from atlas.panelgeneration import VariantPanel


def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]

import argparse
parser = argparse.ArgumentParser(
    description='Parse VCF and upload variants to DB')
parser.add_argument(
    'reference_filepath',
    metavar='reference_filepath',
    type=str,
    help='reference_filepath')
parser.add_argument(
    '--db_name',
    metavar='db_name',
    type=str,
    help='db_name',
    default="tb")
parser.add_argument(
    '--kmer',
    metavar='kmer',
    type=int,
    help='kmer length',
    default=31)
parser.add_argument('--force', default=False, action="store_true")
parser.add_argument('-v', '--verbose', default=False, action="store_true")
args = parser.parse_args()

if args.verbose:
    LOGGER.setLevel(level=logging.DEBUG)
else:
    LOGGER.setLevel(level=logging.ERROR)

connect('atlas-%s' % (args.db_name))


def get_context(pos):
    context = []
    for variant in Variant.objects(
            start__ne=pos,
            start__gt=pos - args.kmer,
            start__lt=pos + args.kmer):
        for split_variant in variant.split():
            context.append(split_variant)
    return context


def called_variant_list_to_panel_generated_variant_list(variants_in):
    variants_out = []
    for var in variants_in:
        for alt in var.alternate_bases:
            variants_out.append(Variant(var.reference_bases, var.start, alt))
    return variants_out


def seen_together(variants):
    # Takes a list of variants.
    # Returns a list of variants that appear together (in the same variant set)
    variant_to_samples = {}
    for variant in variants:
        variant_to_samples[variant] = variant.seen_in_samples()

    samples_counter = Counter(flatten(variant_to_samples.values()))
    samples_seen_more_than_once = [
        k for k,
        v in samples_counter.iteritems() if v > 1]  # 2 as we also have global var
    contexts = []
    for sample in samples_seen_more_than_once:
        vars_together = []
        for variant, samples in variant_to_samples.items():
            if sample in samples:
                vars_together.append(variant)
        if vars_together not in contexts:
            contexts.append(vars_together)
            variants = [var for var in variants if var not in vars_together]
    for var in variants:
        contexts.append([var])
    return contexts + [[]]


def make_variant_probe(variant):
    context = get_context(variant.start)
    variant_probe = None
    contexts_seen_together = seen_together(context)
    alts = []
    build_success = True
    for context in contexts_seen_together:
        if len(context) <= 5:
            try:
                panel = al.create(variant, context)
                ref = panel.ref
                panel.alts
                if variant_probe is not None:
                    variant_probe.alts.extend(panel.alts)
                else:
                    variant_probe = panel
            except ValueError as e:
                print e
                build_success = False
    variant_probe.alts = unique(variant_probe.alts)
    return variant_probe


al = AlleleGenerator(
    reference_filepath=args.reference_filepath,
    kmer=args.kmer)
for variant in Variant.snps():
    variant_panel = make_variant_probe(variant)
    sys.stdout.write(
        ">ref-%s?num_alts=%i&ref=%s\n" %
        (variant_panel.variant.var_hash, len(
            variant_panel.alts), variant_panel.variant.reference.id))
    sys.stdout.write("%s\n" % variant_panel.ref)
    for a in variant_panel.alts:
        sys.stdout.write(">alt-%s\n" % variant_panel.variant.var_hash)
        sys.stdout.write("%s\n" % a)
