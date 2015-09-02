#! /usr/bin/env python
import time
import json
import pickle
from treeplacing import Node
from treeplacing import Leaf
from treeplacing import Placer
from mongoengine import connect
import argparse
parser = argparse.ArgumentParser(description='Genotype a sample based on kmer coverage on alleles')
parser.add_argument('sample', metavar='sample', type=str, help='sample id')
parser.add_argument('db_name', metavar='db_name', type=str, help='db_name', default="default")
parser.add_argument('kmer', metavar='kmer', type=int, help='kmer size')
args = parser.parse_args()

connect('atlas-%s-%i' % (args.db_name ,args.kmer))

start = time.clock()

with open ("der_tree.json", 'r') as infile:
    tree = json.load(infile)

# def walk(node):
#     children = []
#     for child in node['children']:
#         if child["type"] == "node":
#             childrens = walk(child)
#             children.append(Node(children = childrens))
#         else:
#             children.append(Leaf(child["name"]))
#     return children

# root = Node(children = walk(tree))

# pickle.dump( root, open( "root.p", "wb" ) ) 
root = pickle.load( open( "root.p", "rb" ) )
print root
neighbours = Placer(root).place(args.sample)
if type(neighbours) is list:
	print "Nearest Neighbours are %s" % ",".join(neighbours)
else:
	print "Nearest Neighbour is %s" % neighbours
print "Time taken to search %s seconds " % str((time.clock() - start))