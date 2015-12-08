#! /usr/bin/env python
from __future__ import print_function

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import csv

from mongoengine import connect

from collections import Counter

from Bio.Data import CodonTable 

from atlas.panelgeneration import AlleleGenerator
from atlas.panelgeneration import Variant

from atlas.vcf2db import split_var_name
from atlas.vcf2db import VariantFreq
from atlas.vcf2db import Variant as CalledVariant

from atlas.genes import GeneAminoAcidChangeToDNAVariants

import argparse

parser = argparse.ArgumentParser(description='Parse VCF and upload variants to DB')
parser.add_argument('reference_filepath', metavar='reference_filepath', type=str, help='reference_filepath')
parser.add_argument('-v','--variant', type=str, action='append', help='Variant in DNA positions e.g. A1234T', default = [])
parser.add_argument('-f','--file', type=str, help='File containing variants as rows A1234T')
parser.add_argument('-g','--genbank', type=str, help='Genbank file containing genes as features')
parser.add_argument('--db_name', metavar='db_name', type=str, help='db_name', default="tb")
parser.add_argument('--kmer', metavar='kmer', type=int, help='kmer length', default = 31)
# parser.add_argument('--alphabet', metavar='alphabet', type=str, help='DNA or PROT variants', choices = ["DNA", "PROT"])
parser.add_argument('-q', '--quiet', default = False, action = "store_true")
parser.add_argument('--mykrobe', default = False, action = "store_true")
args = parser.parse_args()

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

mutations = []
## Check if variants are in aminoacid space
class Mutation(object):

    def __init__(self, dna_var, gene = None, mut = None):
        self.dna_var = dna_var
        self.gene = gene
        tmp, self.location, tmp= split_var_name(mut)
        self.ref, tmp, self.alt = split_var_name(dna_var)

    @property
    def mut(self):
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        r = standard_table.forward_table.get(self.ref, self.ref)
        a = standard_table.forward_table.get(self.alt, self.alt)
        return "".join([r, str(self.location), a])

if args.genbank:
    aa2dna = GeneAminoAcidChangeToDNAVariants(args.reference_filepath, args.genbank)
    if args.file:
        with open(args.file, 'r') as infile:
            reader = csv.reader(infile, delimiter = "\t")
            for row in reader:
                gene, mutation = row
                for dna_var in aa2dna.get_variant_names(gene, mutation):
                    mutations.append(Mutation(dna_var, gene = gene, mut = mutation))
    else:
        for variant in args.variant:
            gene, mutation = variant.split("_")
            for dna_var in aa2dna.get_variant_names(gene, mutation):
                mutations.append(Mutation(dna_var, gene = gene, mut = mutation))
else:    
    if args.file:
        with open(args.file, 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                mutations.append(Mutation(dna_var = row[0]))
    else:
        mutations.extend(Mutation(v) for v in args.variants)

def get_context(pos):
    context = []
    for vft in VariantFreq.objects(start__ne = pos, start__gt = pos - args.kmer, start__lt = pos + args.kmer):
        for alt in vft.alternate_bases:
            context.append( Variant(vft.reference_bases, vft.start , alt) )
    return context

def seen_together(variants):
    ## Takes a list of variants. 
    ## Returns a list of variants that appear together (in the same variant set)
    ## TODO - would be better if seen in same sample (accross several variant sets)
    variant_name_hashes = [v.name_hash for v in variants]
    variant_set_counter = Counter([v.variant_set.id for v in CalledVariant.objects(name_hash__in = variant_name_hashes)])
    variant_sets = [k for k,v in variant_set_counter.iteritems() if v > 1]
    contexts = []
    for vs in variant_sets:
        vars_together = [v for v in CalledVariant.objects(name_hash__in = variant_name_hashes, variant_set = vs)]
        if not vars_together in contexts:
            contexts.append(called_variant_list_to_panel_generated_variant_list(vars_together))
            variants = [var for var in variants if var not in vars_together] 
    for var in variants:
        contexts.append([var])
    if not contexts:
        return [[]]
    else:
        return contexts

def make_panels(var):
    reference_bases, start, alt = split_var_name(var)
    alternate_bases = alt.split('/')
    context = get_context(start)
    panels = []
    contexts_seen_together = seen_together(context)    
    for alt in alternate_bases:
        variant = Variant(reference_bases, start , alt)
        for context in contexts_seen_together:
            if len(context) <= 5:
                try:
                    panel = al.create(variant, context)
                except ValueError, e:
                    sys.stderr.write(str(e))
                    panel = al.create(variant)
                vo = "".join([reference_bases, str(start), alt])
                panels.append((vo, panel))
    return panels

al = AlleleGenerator(reference_filepath = args.reference_filepath, kmer = args.kmer)
if args.mykrobe:
    for mut in mutations:
        var = mut.dna_var
        panels = make_panels(var)
        for name, variant_panel in panels:
            sys.stdout.write(">ref_%s_sub-%i-%s\n" % (mut.mut, len(variant_panel.alts), mut.gene))
            sys.stdout.write("%s\n" % variant_panel.ref)
            for i,a in enumerate(variant_panel.alts):
                sys.stdout.write(">alt_%s_alt-%i-%s\n" % (mut.mut, i + 1, mut.gene))
                sys.stdout.write("%s\n" % a)  
else:
    for mut in mutations:
        var = mut.dna_var
        panels = make_panels(var)
        for name, variant_panel in panels:
            if mut.gene:
                sys.stdout.write(">ref-%s?num_alts=%i&gene=%s&mut=%s\n" % (name, len(variant_panel.alts), mut.gene, mut.mut ))
            else:
                sys.stdout.write(">ref-%s?num_alts=%i\n" % (name, len(variant_panel.alts)))
            sys.stdout.write("%s\n" % variant_panel.ref)
            for a in variant_panel.alts:
                sys.stdout.write(">alt-%s\n" % name)
                sys.stdout.write("%s\n" % a)        
