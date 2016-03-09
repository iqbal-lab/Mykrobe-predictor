from __future__ import print_function
import csv
import os
import sys
import logging
from mongoengine import connect
from mongoengine.connection import ConnectionError
from pymongo.errors import ServerSelectionTimeoutError
from Bio.Seq import Seq

from atlas.panelgeneration import AlleleGenerator
from atlas.panelgeneration import make_variant_probe
from atlas.schema import Variant
from atlas.utils import split_var_name
from atlas.annotation.genes import GeneAminoAcidChangeToDNAVariants
from atlas.panelgeneration.models import Mutation

# logging = logging.getLogger(__name__)
# logging.setLevel(level=logging.DEBUG)


def run(parser, args):
    DB = connect('atlas-%s' % (args.db_name))
    if DB is not None:
        try:
            Variant.objects()
        except (ServerSelectionTimeoutError, ConnectionError):
            DB = None
            logging.warning(
                "Could not connect to database. Continuing without using genetic backgrounds")
    mutations = []
    reference = os.path.basename(args.reference_filepath).split('.fa')[0]
    if args.genbank:
        aa2dna = GeneAminoAcidChangeToDNAVariants(
            args.reference_filepath,
            args.genbank)
        if args.file:
            with open(args.file, 'r') as infile:
                reader = csv.reader(infile, delimiter="\t")
                for row in reader:
                    gene, mutation, alphabet = row
                    if alphabet == "DNA":
                        protein_coding_var = False
                    else:
                        protein_coding_var = True
                    for var_name in aa2dna.get_variant_names(
                            gene, mutation, protein_coding_var):
                        mutations.append(
                            Mutation(reference=reference,
                                     var_name=var_name,
                                     gene=aa2dna.get_gene(gene),
                                     mut=mutation))
        else:
            for variant in args.variant:
                gene, mutation = variant.split("_")
                for var_name in aa2dna.get_variant_names(gene, mutation):
                    mutations.append(
                        Mutation(reference=reference,
                                 var_name=var_name,
                                 gene=gene,
                                 mut=mutation))
    else:
        if args.file:
            with open(args.file, 'r') as infile:
                reader = csv.reader(infile)
                for row in reader:
                    mutations.append(
                        Mutation(
                            reference=reference,
                            var_name=row[0]))
        else:
            mutations.extend(Mutation(reference=reference, var_name=v)
                             for v in args.variants)

    al = AlleleGenerator(
        reference_filepath=args.reference_filepath,
        kmer=args.kmer)
    for mut in mutations:
        variant_panel = make_variant_probe(al, mut.variant, args.kmer, DB=DB)
        if variant_panel is not None:
            if mut.gene:
                sys.stdout.write(
                    ">ref-%s?num_alts=%i&gene=%s&mut=%s&ref=%s\n" %
                    (mut.variant.var_name, len(
                        variant_panel.alts), mut.gene.name, mut.mut, mut.reference))
            else:
                sys.stdout.write(
                    ">ref-%s?num_alts=%i\n" %
                    (mut.variant.var_name, len(
                        variant_panel.alts)))
            sys.stdout.write("%s\n" % variant_panel.ref)
            for a in variant_panel.alts:
                sys.stdout.write(">alt-%s\n" % mut.mut)
                sys.stdout.write("%s\n" % a)
        else:
            logging.warning(
                "All variants failed for %s_%s - %s" %
                (mut.gene, mut.mut, mut.variant))
