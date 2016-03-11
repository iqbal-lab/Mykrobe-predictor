# Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.typing import CoverageParser
from atlas.metagenomics import AMRSpeciesPredictor
from atlas.version import __version__

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
        tmp_dir=args.tmp,
        skeleton_dir=args.skeleton_dir,
        mccortex31_path=args.mccortex31_path)
    cp.run()

    if args.expected_depth is None:
        args.expected_depth = cp.estimate_depth()

    base_json = {args.sample: {}}
    base_json[args.sample]["panels"] = args.panels
    base_json[args.sample]["files"] = args.seq
    base_json[args.sample]["kmer"] = args.kmer
    base_json[args.sample]["version"] = __version__
    gt = Genotyper(sample=args.sample, expected_depths=[args.expected_depth],
                   variant_covgs=cp.variant_covgs,
                   gene_presence_covgs=cp.covgs["presence"],
                   base_json=base_json,
                   contamination_depths=[],
                   force_gt=args.force_gt)
    gt.run()
    cp.remove_temporary_files()
    print(json.dumps(gt.out_json, indent=4))
