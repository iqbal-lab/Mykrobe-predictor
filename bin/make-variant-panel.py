#! /usr/bin/env python
from __future__ import print_function

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
import csv
from mongoengine import connect
from Bio.Data import CodonTable
from Bio.Seq import Seq

from atlas.schema import Variant
from atlas.panelgeneration import AlleleGenerator
from atlas.panelgeneration import make_variant_probe
from atlas.utils import split_var_name
from atlas.annotation.genes import GeneAminoAcidChangeToDNAVariants

import argparse

parser = argparse.ArgumentParser(
    description='Parse VCF and upload variants to DB')
parser.add_argument(
    'reference_filepath',
    metavar='reference_filepath',
    type=str,
    help='reference_filepath')
parser.add_argument(
    '-v',
    '--variant',
    type=str,
    action='append',
    help='Variant in DNA positions e.g. A1234T',
    default=[])
parser.add_argument(
    '-f',
    '--file',
    type=str,
    help='File containing variants as rows A1234T')
parser.add_argument(
    '-g',
    '--genbank',
    type=str,
    help='Genbank file containing genes as features')
parser.add_argument(
    '--db_name',
    metavar='db_name',
    type=str,
    help='db_name',
    default="tb")
parser.add_argument('-k', '--kmer', type=int, help='kmer length', default=31)
# parser.add_argument('--alphabet', metavar='alphabet', type=str, help='DNA or PROT variants', choices = ["DNA", "PROT"])
parser.add_argument('-q', '--quiet', default=False, action="store_true")
parser.add_argument('--mykrobe', default=False, action="store_true")
args = parser.parse_args()

DB_NAME = 'atlas-%s' % (args.db_name)
connect(DB_NAME)

mutations = []
# Check if variants are in aminoacid space


class Mutation(object):

    def __init__(self, var_name, gene=None, mut=None, reference = os.path.basename(args.reference_filepath).split('.')[0]):
        self.var_name = var_name
        self.gene = gene
        if mut:
            tmp, self.start, tmp = split_var_name(mut)
        self.ref, tmp, self.alt = split_var_name(var_name)
        self.standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.reference = reference

    @property
    def mut(self):
        if self.gene.forward:
            ref = self.ref
            alt = self.alt
        else:
            ref = str(Seq(self.ref).reverse_complement())
            alt = str(Seq(self.alt).reverse_complement())
        r = self.standard_table.forward_table.get(ref, ref)
        a = self.standard_table.forward_table.get(alt, alt)
        return "".join([r, str(self.start), a])

    @property
    def variant(self):
        ref, start, alt = split_var_name(self.var_name)
        return Variant.create(variant_sets=None, start=int(start),
                            end = 0, reference_bases=ref,
                            alternate_bases=[alt],
                            reference=self.reference)
    

def run(parser, args):
    if args.genbank:
        aa2dna = GeneAminoAcidChangeToDNAVariants(
            args.reference_filepath,
            args.genbank)
        if args.file:
            with open(args.file, 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                for row in reader:
                    gene, mutation = row
                    for var_name in aa2dna.get_variant_names(gene, mutation):
                        mutations.append(
                            Mutation(
                                var_name,
                                gene=aa2dna.get_gene(gene),
                                mut=mutation))
        else:
            for variant in args.variant:
                gene, mutation = variant.split("_")
                for var_name in aa2dna.get_variant_names(gene, mutation):
                    mutations.append(Mutation(var_name, gene=gene, mut=mutation))
    else:
        if args.file:
            with open(args.file, 'r') as infile:
                reader = csv.reader(infile)
                for row in reader:
                    mutations.append(Mutation(var_name=row[0]))
        else:
            mutations.extend(Mutation(v) for v in args.variants)

    al = AlleleGenerator(
        reference_filepath=args.reference_filepath,
        kmer=args.kmer)
    if args.mykrobe:
        for mut in mutations:
            variant = mut.variant
            panels = make_variant_probe(al, variant, args.kmer)
            for name, variant_panel in panels:
                sys.stdout.write(
                    ">ref_%s_sub-%i-%s\n" %
                    (mut.mut, len(
                        variant_panel.alts), mut.gene.name))
                sys.stdout.write("%s\n" % variant_panel.ref)
                for i, a in enumerate(variant_panel.alts):
                    sys.stdout.write(
                        ">alt_%s_alt-%i-%s\n" %
                        (mut.mut, i + 1, mut.gene))
                    sys.stdout.write("%s\n" % a)
    else:
        for mut in mutations:
            variant_panel = make_variant_probe(al, mut.variant, args.kmer)
            # for name, variant_panel in panels:
            if mut.gene:
                sys.stdout.write(
                    ">ref-%s?num_alts=%i&gene=%s&mut=%s&ref=%s\n" %
                    (mut.variant.var_name, len(
                        variant_panel.alts), mut.gene.name, mut.mut, os.path.basename(
                        args.reference_filepath).split('.')[0]))
            else:
                sys.stdout.write(
                    ">ref-%s?num_alts=%i\n" %
                    (mut.variant.var_name, len(
                        variant_panel.alts)))
            sys.stdout.write("%s\n" % variant_panel.ref)
            for a in variant_panel.alts:
                sys.stdout.write(">alt-%s\n" % mut.mut)
                sys.stdout.write("%s\n" % a)
run(parser, args)
