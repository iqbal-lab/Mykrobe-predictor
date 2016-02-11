import os
import hashlib
import re
import json

def check_args(args):
    if args.db_name is None:
        args.db_name = os.environ.get("DB_NAME")
    if args.db_name is None:
        raise ValueError("db_name needs to be set. Either run with --db_name :db_name or export DB_NAME=:db_name")
    if args.kmer is None:
        args.kmer = os.environ.get("KMER_SIZE")
    if args.kmer is None:
        raise ValueError("kmer needs to be set. Either run with --kmer :kmer_size or export KMER_SIZE=:kmer_size")
    else:
        args.kmer = int(args.kmer)  
    return args

def make_hash(s):
    return hashlib.sha256(s.encode("ascii", errors="ignore")).hexdigest()   

def split_var_name(name):
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/]+)", name, re.I).groups()
    return items[0],int(items[1]),items[2] 

def unique(l):
    seen = set()
    return [x for x in l if x not in seen and not seen.add(x)]

def flatten(l):
    return [item for sublist in l for item in sublist]

def get_params(url):
    params = {}
    try:
        p_str = url.split("?")[1]
    except IndexError:
        return params
    p_str = p_str.split(" ")[0]
    p_str = p_str.split('&')
    for p in p_str:
        k,v = p.split("=")
        params[k] = v
    return params

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def load_json(f):
    with open(f, 'r') as infile:
        return json.load(infile)        