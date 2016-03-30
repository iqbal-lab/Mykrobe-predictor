from __future__ import print_function
import logging
from pprint import pprint
import json
import os
from mykatlas.utils import check_args
from mykatlas.typing import CoverageParser
from mykatlas.typing import Genotyper
from mykrobe.predict import TBPredictor
from mykrobe.predict import StaphPredictor
from mykrobe.predict import GramNegPredictor
from mykrobe.metagenomics import AMRSpeciesPredictor
from mykrobe.version import __version__ as predictor_version
from mykatlas.version import __version__ as atlas_version
STAPH_PANELS = ["data/panels/staph-species-160227.fasta.gz",
                "data/panels/staph-amr-bradley_2015.fasta.gz"]

GN_PANELS = [
    "data/panels/gn-amr-genes",
    "data/panels/Escherichia_coli",
    "data/panels/Klebsiella_pneumoniae",
    "data/panels/gn-amr-genes-extended"]


def run(parser, args):
    base_json = {args.sample: {}}
    args = parser.parse_args()
    hierarchy_json_file = None
    if args.panel is not None:
        if args.panel == "bradley-2015":
            TB_PANELS = [
                "data/panels/tb-species-160330.fasta.gz",
                "data/panels/tb-amr-bradley_2015.fasta.gz"]
        elif args.panel == "walker-2015":
            TB_PANELS = [
                "data/panels/tb-species-160330.fasta.gz",
                "data/panels/tb-amr-walker_2015.fasta.gz"]

    if not args.species:
        panels = TB_PANELS + GN_PANELS + STAPH_PANELS
        panel_name = "tb-gn-staph-amr"

    elif args.species == "staph":
        panels = STAPH_PANELS
        panel_name = "staph-amr"
        # hierarchy_json_file = "data/phylo/saureus_hierarchy.json"

    elif args.species == "tb":
        panels = TB_PANELS
        panel_name = "tb-amr"
        hierarchy_json_file = "data/phylo/mtbc_hierarchy.json"
    elif args.species == "gn":
        panels = GN_PANELS
        panel_name = "gn-amr"
    logging.info("Running AMR prediction with panels %s" % ", ".join(panels))
    base_json[args.sample]["panels"] = panels
    base_json[args.sample]["files"] = args.seq
    base_json[args.sample]["kmer"] = args.kmer
    base_json[args.sample]["version"] = {}
    base_json[args.sample]["version"]["mykrobe-predictor"] = predictor_version
    base_json[args.sample]["version"]["mykrobe-atlas"] = atlas_version
    # Get real paths for panels
    panels = [
        os.path.realpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                f)) for f in panels]
    if hierarchy_json_file is not None:
        hierarchy_json_file = os.path.realpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                hierarchy_json_file))
    # Run Cortex
    cp = CoverageParser(
        sample=args.sample,
        panel_file_paths=panels,
        seq=args.seq,
        kmer=args.kmer,
        force=args.force,
        verbose=False,
        tmp_dir=args.tmp,
        skeleton_dir=args.skeleton_dir,
        mccortex31_path=args.mccortex31_path)
    cp.run()
    # Detect species
    species_predictor = AMRSpeciesPredictor(
        phylo_group_covgs=cp.covgs.get(
            "complex",
            cp.covgs.get(
                "phylo_group",
                {})),
        sub_complex_covgs=cp.covgs.get(
            "sub-complex",
            {}),
        species_covgs=cp.covgs["species"],
        lineage_covgs=cp.covgs.get(
            "sub-species",
            {}),
        base_json=base_json[args.sample],
        hierarchy_json_file=hierarchy_json_file)
    species_predictor.run()

    # ## AMR prediction

    depths = []
    Predictor = None
    if species_predictor.is_saureus_present():
        depths = [species_predictor.out_json["phylogenetics"]
                  ["phylo_group"]["Staphaureus"]["median_depth"]]
        Predictor = StaphPredictor
    elif species_predictor.is_mtbc_present():
        depths = [species_predictor.out_json["phylogenetics"]["phylo_group"][
            "Mycobacterium_tuberculosis_complex"]["median_depth"]]
        Predictor = TBPredictor
    elif species_predictor.is_gram_neg_present():
        Predictor = GramNegPredictor
        try:
            depths = [species_predictor.out_json["phylogenetics"][
                "species"]["Klebsiella_pneumoniae"]["median_depth"]]
        except KeyError:
            depths = [species_predictor.out_json["phylogenetics"]
                      ["species"]["Escherichia_coli"]["median_depth"]]
    # pprint (species_predictor.out_json["phylogenetics"]["species"])
    # Genotype
    q = args.quiet
    args.quiet = True
    if depths:
        gt = Genotyper(sample=args.sample, expected_depths=depths,
                       variant_covgs=cp.variant_covgs,
                       gene_presence_covgs=cp.covgs["presence"],
                       base_json=base_json,
                       contamination_depths=[],
                       include_hom_alt_calls=True)
        gt.run()
    args.quiet = q
    if Predictor is not None:
        predictor = Predictor(variant_calls=gt.variant_calls,
                              called_genes=gt.gene_presence_covgs,
                              base_json=base_json[args.sample])
        predictor.run()
    cp.remove_temporary_files()
    print(json.dumps(base_json, indent=4))
