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
    cp = CoverageParser(
        sample=args.sample,
        panel_file_paths=args.panels,
        seq=args.seq,
        kmer=args.kmer,
        force=args.force,
        verbose=verbose,
        skeleton_dir=args.tmp)
    cp.run()
    base_json = {args.sample: {}}
    gt = Genotyper(sample=args.sample, expected_depths=[100],
                   variant_covgs=cp.variant_covgs,
                   gene_presence_covgs=cp.covgs["presence"],
                   base_json=base_json,
                   contamination_depths=[])
    gt.run()
    print(json.dumps(gt.out_json, indent=4))
