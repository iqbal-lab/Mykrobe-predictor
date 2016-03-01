from __future__ import print_function

# Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import CoverageParser
from atlas.typing import Genotyper
from atlas.pheno import TBPredictor
from atlas.pheno import StaphPredictor
from atlas.pheno import GramNegPredictor
from atlas.metagenomics import AMRSpeciesPredictor
from pprint import pprint
import json

STAPH_PANELS = ["Coagneg",
                "Staphaureus",
                "Saureus",
                "Sepidermidis",
                "Shaemolyticus",
                "Sother",
                "staph-amr-genes",
                "staph-amr-mutations"]
GN_PANELS = [
    "gn-amr-genes",
    "Escherichia_coli",
    "Klebsiella_pneumoniae",
    "gn-amr-genes-extended"]
TB_PANELS = [
    "data/panels/tb-species-160227.fasta",
    "data/panels/tb-amr-walker_2015.fasta"]


def run(parser, args):
    base_json = {args.sample: {}}
    args = parser.parse_args()
    if not args.species:
        panels = TB_PANELS + GN_PANELS + STAPH_PANELS
        panel_name = "tb-gn-staph-amr"
    elif args.species == "staph":
        panels = STAPH_PANELS
        panel_name = "staph-amr"
    elif args.species == "tb":
        panels = TB_PANELS
        panel_name = "tb-amr"
    elif args.species == "gn":
        panels = GN_PANELS
        panel_name = "gn-amr"
    # Run Cortex
    cp = CoverageParser(
        sample=args.sample,
        panel_file_paths=panels,
        seq=args.seq,
        kmer=args.kmer,
        force=args.force,
        verbose=False)
    cp.run()
    # print (cp.covgs["species"])
    # Detect species
    species_predictor = AMRSpeciesPredictor(
        phylo_group_covgs=cp.covgs.get(
            "complex",
            {}),
        sub_complex_covgs=cp.covgs.get(
            "sub-complex",
            {}),
        species_covgs=cp.covgs["species"],
        lineage_covgs=cp.covgs.get(
            "sub-species",
            {}),
        base_json=base_json)
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
        predictor = Predictor(typed_variants=gt.variant_calls,
                              called_genes=gt.gene_presence_covgs,
                              base_json=base_json[args.sample])
        predictor.run()

    print(json.dumps(base_json, indent=4))
