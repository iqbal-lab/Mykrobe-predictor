#!/usr/bin/env python
# encoding: utf-8
"""
compare.py
Created on: 2016-03-30
Author: Phelim Bradley

Compares the JSON output of mykrobe predictor

"""

import argparse

parser = argparse.ArgumentParser(description='Compares the JSON output of mykrobe predictor')
parser.add_argument('truth', metavar='truth', type=str, help='truth')
parser.add_argument('analyses', metavar='analyses', type=str, help='analyses', nargs='+')
args = parser.parse_args()

