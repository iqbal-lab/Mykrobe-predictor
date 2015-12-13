## Read the kmer counts into a hash
from atlas.utils import check_args
from atlas.typing import Genotyper

def run(parser, args):
    args = parser.parse_args()
    check_args(args)  
    Genotyper(args).run()
