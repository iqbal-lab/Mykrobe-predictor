#! /usr/bin/env python
import json
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

with open ("der_tree.json", 'r') as infile:
    tree = json.load(infile)


def walk(node):
    children = []
    for child in node['children']:
        if child["type"] == "node":
            childrens = walk(child)
            children.append(Node(children = childrens))
        else:
            children.append(Leaf(child["name"]))
    return children

root = Node(children = walk(tree))
print root
print root.other_child(root.children[0])

print Placer(root).place(args.sample)