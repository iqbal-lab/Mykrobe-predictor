from __future__ import print_function
## Using seed kmers from a gene panel return a DFS through the graph
import sys
sys.path.append('/home/phelimb/git/atlas-core')
from atlas.cortex.server import WebServer
from atlas.cortex.server import McCortexQuery
from atlas.cortex.server import GraphWalker
from atlas.utils import get_params
import socket
import threading
import json
from pprint import pprint
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
from Bio import SeqIO

import argparse
parser = argparse.ArgumentParser(description='Add length argument to panel')
parser.add_argument('dna_fasta', metavar='dna_fasta', type=str, help='dna_fasta')
parser.add_argument('--prot_fasta', metavar='prot_fasta', type=str, help='prot_fasta')
parser.add_argument('-f', '--ctx', metavar='ctx', type=str, help='cortex graph binary')
parser.add_argument('-k','--kmer_size', metavar='kmer_size', type=int,
                   help='kmer_size', default = 31)
parser.add_argument('-p','--port', metavar='port', type=int,
                   help='port', default = None)
args = parser.parse_args()



def get_open_port():
	if args.port:
		return args.port
	else:
		s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		s.bind(("",0))
		s.listen(1)
		port = s.getsockname()[1]
		s.close()
		return port

class PathDetails(object):

	def __init__(self, start_kmer, last_kmer, length, v = ""):
		self.start_kmer = start_kmer
		self.last_kmer = last_kmer
		self.length = length
		self.version = v

	def __eq__(self, pd):
		return self.start_kmer == pd.start_kmer and self.last_kmer == pd.last_kmer and self.length == self.length

genes = {}
skip_list = {"tem"  : ["191", "192"],
			 "oxa"  :["12", "14", "33"],
			 "shv" : ["12", "6"]
			}

def get_repeat_kmers(record, k):
	## Process repeat kmers
	kmers = {}
	kmers_seq = []
	# repeat_kmers = {}
	for i in range(len(record.seq) - k + 1):
		_kmers = [str(record.seq[i:i+k]), str(record.seq[i:i+k].reverse_complement())]
		for kmer in _kmers:
			kmers_seq.append(kmers_seq)
			if kmers.has_key(kmer):
			  c = max(kmers[kmer].keys()) + 1
			  kmers[kmer][c] = str(record.seq[i + 1 :i + k + 1])
			else:
			  kmers[kmer] = {}
			  kmers[kmer][1] = str(record.seq[i + 1 :i + k + 1])
	# repeat_kmers = {}
	# for kmer, count in kmers.items():
	# 	if len(count.keys()) > 1:
	# 		repeat_kmers[kmer] = count
	return kmers



with open(args.dna_fasta, 'r') as infile:
	for i, record in enumerate(SeqIO.parse(infile, "fasta")):
		repeat_kmers = get_repeat_kmers(record, args.kmer_size)
		params = get_params(record.id)
		gene_name = params.get("name", i)
		version = params.get("version", i)
		start_kmer = str(record.seq)[:args.kmer_size]
		last_kmer = str(record.seq)[-args.kmer_size:]
		if not version in skip_list.get(gene_name, []):
			if not gene_name in genes:
				genes[gene_name] = {}
				genes[gene_name]["pathdetails"] = []
				genes[gene_name]["known_kmers"] = ""
			pd = PathDetails(start_kmer, last_kmer, len(record.seq), v = version)
			# print (gene_name, version, len(record.seq))
			if pd not in genes[gene_name]["pathdetails"]:
				genes[gene_name]["pathdetails"].append(pd)
			else:
				j = genes[gene_name]["pathdetails"].index(pd)
				genes[gene_name]["pathdetails"][j].version += "-%s" % version
		genes[gene_name]["repeat_kmers"] = repeat_kmers
		genes[gene_name]["known_kmers"] += "%sN" % str(record.seq)

with open("gene.tmp.json", "w") as outf:
	json.dump(repeat_kmers, outf, sort_keys = False, indent = 4)

wb = None
if args.port is None:
	if args.cts is None:
		raise ValueError("Require either port or binary")
		port =  (get_open_port())
		logger.debug("Running server on port %i " % port)
		wb = WebServer(port, args = [ "-q",  args.ctx ] )
		logger.debug("Loading binary")
		wb.start()
		## Serve on a thread
		logger.debug("Starting sever")
		server = threading.Thread(target=wb.serve)
		server.start()
else:
	port = args.port


logger.debug("Walking the graph")
out_dict = {}
gw = GraphWalker(port = port, kmer_size = args.kmer_size)

def get_paths_for_gene(gene_name, gene_dict, gw):
	paths = {}
	for pd in gene_dict["pathdetails"]:
		p = gw.breath_first_search(N = pd.length, seed = pd.start_kmer,
		                            end_kmers = [pd.last_kmer],
		                            known_kmers = gene_dict["known_kmers"], 
		                            repeat_kmers = gene_dict["repeat_kmers"])
		print (p)
		if p:
			paths[pd.version] = p
	return paths
for gene_name, gene_dict in genes.items():
	paths = get_paths_for_gene(gene_name, gene_dict, gw)
	out_dict[gene_name] = paths
print (json.dumps(out_dict, sort_keys = False, indent = 4))
logger.info("Cleaning up")
if wb is not None:
	wb.stop()




