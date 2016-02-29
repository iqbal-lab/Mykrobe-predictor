#! /usr/bin/env python

# import sys
# import os
# sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")

import argparse


args = parser.parse_args()

DB_NAME = 'atlas-%s' % (args.db_name)
connect(DB_NAME)

mutations = []
# Check if variants are in aminoacid space
