## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper
from atlas.typing import CoverageParser
import json
def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    cp = CoverageParser(args, panels = ["gn-plasmids"], verbose = False)
    cp.run()    
    gt = Genotyper(args, depths = [100],
                variant_covgs = cp.covgs["variant"],
                gene_presence_covgs = cp.covgs["presence"],
                verbose = False, 
                base_json = cp.out_json,
                contamination_depths = [])
    gt.run()
    print(json.dumps(cp.out_json, indent = 4))    