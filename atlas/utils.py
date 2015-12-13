import os
import hashlib
import re

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