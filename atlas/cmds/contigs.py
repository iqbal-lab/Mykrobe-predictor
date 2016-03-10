from __future__ import print_function
import sys
from atlas.cortex.server import WebServer
from atlas.cortex.server import McCortexQuery
from atlas.cortex.server import GraphWalker
from atlas.cortex.server import query_mccortex
from atlas.utils import get_params
import socket
import json
from pprint import pprint
import logging
from Bio import SeqIO
import argparse

sys.path.append('/home/phelimb/git/atlas-core')

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Using seed kmers from a gene panel return a DFS through the graph

parser = argparse.ArgumentParser(description='Add length argument to panel')
parser.add_argument(
    'dna_fasta',
    metavar='dna_fasta',
    type=str,
    help='dna_fasta')
parser.add_argument(
    '--prot_fasta',
    metavar='prot_fasta',
    type=str,
    help='prot_fasta')
parser.add_argument(
    '-f',
    '--ctx',
    metavar='ctx',
    type=str,
    help='cortex graph binary')
parser.add_argument('-k', '--kmer_size', metavar='kmer_size', type=int,
                    help='kmer_size', default=31)
parser.add_argument('-p', '--port', metavar='port', type=int,
                    help='port', default=None)
args = parser.parse_args()


# def get_open_port():
# 	if args.port:
# 		return args.port
# 	else:
# 		s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
# 		s.bind(("",0))
# 		s.listen(1)
# 		port = s.getsockname()[1]
# 		s.close()
# 		return port

class PathDetails(object):

    def __init__(self, start_kmer, last_kmer, length, skipped=0, v=""):
        self.start_kmer = start_kmer
        self.last_kmer = last_kmer
        self.length = length
        self.skipped = skipped
        self.version = v
        self.repeat_kmers = {}

    def __eq__(self, pd):
        return self.start_kmer == pd.start_kmer and self.last_kmer == pd.last_kmer and self.length == pd.length and self.skipped == pd.skipped

    def set_repeat_kmers(self, repeat_kmers):
        if not self.repeat_kmers:
            self.repeat_kmers = repeat_kmers
        else:
            raise ValueError("Already set repeat kmers")

genes = {}
skip_list = {"tem": ["191", "192"],
             "oxa": ["12", "14", "33"],
             "shv": ["12", "6"]
             }


def get_repeat_kmers(record, k):
    # Process repeat kmers
    kmers = {}
    kmers_seq = []
    # repeat_kmers = {}
    for i in range(len(record.seq) - k + 1):
        _kmers = [str(record.seq[i:i + k]),
                  str(record.seq[i:i + k].reverse_complement())]
        for kmer in _kmers:
            kmers_seq.append(kmers_seq)
            if kmer in kmers:
                c = max(kmers[kmer].keys()) + 1
                kmers[kmer][c] = str(record.seq[i + 1:i + k + 1])
            else:
                kmers[kmer] = {}
                kmers[kmer][1] = str(record.seq[i + 1:i + k + 1])
    # repeat_kmers = {}
    # for kmer, count in kmers.items():
    # 	if len(count.keys()) > 1:
    # 		repeat_kmers[kmer] = count
    return kmers


def find_start_kmer(seq, mcq, k=args.kmer_size):
    skipped = 0
    for i in range(len(record.seq) - k + 1):
        kmer = seq[i:i + k]
        q = mcq.query(kmer)
        if q.data and q.depth > 1:
            return kmer, skipped
        skipped += 1
    return None, -1


def choose_best_assembly(paths):
    paths.sort(key=lambda x: x["min_non_zero_depth"], reverse=True)
    current_best = paths[0]
    for path in paths[1:]:
        if path["min_non_zero_depth"] < current_best["min_non_zero_depth"]:
            return current_best
        elif path["min_non_zero_depth"] == current_best["min_non_zero_depth"]:
            if path["median_depth"] > current_best["median_depth"]:
                current_best = path
    return current_best


def get_paths_for_gene(gene_name, gene_dict, gw):
    paths = {}
    current_prots = []
    for pd in gene_dict["pathdetails"]:
        path_for_pd = gw.breath_first_search(
            N=pd.length,
            seed=pd.start_kmer,
            end_kmers=[
                pd.last_kmer],
            known_kmers=gene_dict["known_kmers"],
            repeat_kmers=pd.repeat_kmers,
            N_left=pd.skipped)
        # print (p)
        if path_for_pd:
            keep_p = []
            for p in path_for_pd:
                if p["prot"] not in current_prots:
                    keep_p.append(p)
                    current_prots.append(p["prot"])
            if keep_p:
                if len(keep_p) > 1:
                    raise NotImplementedError()
                else:
                    paths[pd.version] = keep_p[0]
    return paths


wb = None
if args.port is None:
    if args.ctx is None:
        raise ValueError(
            "Require either port or binary. e.g. atlas contigs panel.fasta -f sample.ctx")
    else:
        # port =  (get_open_port())
        # logger.debug("Running server on port %i " % port)
        wb = WebServer(port=0, args=[args.ctx])
        logger.debug("Loading binary")
        wb.start()
        # Serve on a thread
        # logger.debug("Starting sever")
        # server = multiprocessing.Process(target=wb.serve, args = (1000,))
        # server.start()
else:
    port = args.port

logger.debug("Walking the graph")
out_dict = {}
gw = GraphWalker(proc=wb.mccortex, kmer_size=args.kmer_size, print_depths=True)


with open(args.dna_fasta, 'r') as infile:
    for i, record in enumerate(SeqIO.parse(infile, "fasta")):
        repeat_kmers = get_repeat_kmers(record, args.kmer_size)
        params = get_params(record.id)
        gene_name = params.get("name", i)
        version = params.get("version", i)
        # print (params, gene_name, version)

        # start_kmer = str(record.seq)[:args.kmer_size]
        last_kmer = str(record.seq)[-args.kmer_size:]
        start_kmer, skipped = find_start_kmer(
            str(record.seq), gw.mcq, args.kmer_size)
        if gene_name not in genes:
            genes[gene_name] = {}
            genes[gene_name]["pathdetails"] = []
            genes[gene_name]["known_kmers"] = ""
        if version not in skip_list.get(gene_name, []) and start_kmer:
            pd = PathDetails(start_kmer, last_kmer, len(record.seq),
                             skipped=skipped, v=version)
            pd.set_repeat_kmers(repeat_kmers)
            # print (gene_name, version, len(record.seq))
            # if pd not in genes[gene_name]["pathdetails"]:
            genes[gene_name]["pathdetails"].append(pd)
            # else:
            # j = genes[gene_name]["pathdetails"].index(pd)
            # _pd = genes[gene_name]["pathdetails"][j]
            # _pd.version += "-%s" % version
            # _pd.set_repeat_kmers(repeat_kmers)

        if gene_name in genes:
            genes[gene_name]["known_kmers"] += "%sN" % str(record.seq)


# logger.debug(len(genes.get("blaZ",{}).get("pathdetails")))
# with open("gene.tmp.json", "w") as outf:
# 	json.dump(genes, outf, sort_keys = False, indent = 4)

for gene_name, gene_dict in genes.items():
    paths = get_paths_for_gene(gene_name, gene_dict, gw)
    if len(paths.keys()) > 1:
        # choose best version
        best_path = choose_best_assembly(paths.values())
    elif len(paths.keys()) == 1:
        best_path = paths.values()[0]
    else:
        best_path = {}
    out_dict[gene_name] = best_path
print (json.dumps(out_dict, sort_keys=False, indent=4))
logger.info("Cleaning up")
if wb is not None:
    wb.stop()
