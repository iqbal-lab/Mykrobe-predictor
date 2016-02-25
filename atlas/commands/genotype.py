# Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.typing import CoverageParser
from atlas.metagenomics import AMRSpeciesPredictor

from pprint import pprint
import json


def run(parser, args):
    args = parser.parse_args()
    # check_args(args)

    # panels = ["tb-species-extended"]

    verbose = True
    cp = CoverageParser(args, panel_file_paths=args.panels, verbose=verbose,
                skeleton_dir = args.tmp)
    cp.run()
    # species_predictor = AMRSpeciesPredictor(
    #     phylo_group_covgs=cp.covgs.get(
    #         "complex",
    #         {}),
    #     sub_complex_covgs=cp.covgs.get(
    #         "sub-complex",
    #         {}),
    #     species_covgs=cp.covgs["species"],
    #     lineage_covgs=cp.covgs.get(
    #         "sub-species",
    #         {}),
    #     base_json=cp.out_json[
    #         args.sample],
    #     verbose=False)
    # species_predictor.run()
    # pprint (species_predictor.out_json["phylogenetics"]["species"])
    base_json = {args.sample : {}}
    gt = Genotyper(args, expected_depths = [100],
                variant_covgs = cp.variant_covgs,
                gene_presence_covgs = cp.covgs["presence"],
                verbose = verbose,
                base_json = base_json,
                contamination_depths = [])
    gt.run()
    # print(json.dumps(cp.out_json, indent=4))
