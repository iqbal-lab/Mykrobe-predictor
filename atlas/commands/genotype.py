## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.typing import CoverageParser
from atlas.metagenomics import AMRSpeciesPredictor

from pprint import pprint
import json

def run(parser, args):
    args = parser.parse_args()
    check_args(args)  

    panels = ["tb-species-extended"]

    verbose = True
    cp = CoverageParser(args, panels = panels, verbose = verbose)
    cp.run() 
    species_predictor = AMRSpeciesPredictor(phylo_group_covgs = cp.covgs.get("complex", {}),
                                            sub_complex_covgs = cp.covgs.get("sub-complex",{}),
                                            species_covgs = cp.covgs["species"],
                                            lineage_covgs = cp.covgs.get("sub-species", {}),
                                            base_json = cp.out_json[args.sample])
    species_predictor.run()
    # pprint (species_predictor.out_json["phylogenetics"]["species"])       
    # gt = Genotyper(args, depths = [100],
    #             variant_covgs = cp.covgs["variant"],
    #             gene_presence_covgs = cp.covgs["presence"],
    #             verbose = verbose, 
    #             base_json = cp.out_json,
    #             contamination_depths = [])
    # gt.run()
    print(json.dumps(species_predictor.out_json["phylogenetics"], indent = 4))    