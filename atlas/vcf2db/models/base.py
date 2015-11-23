import hashlib
import re
def make_hash(s):
	return hashlib.sha256(s.encode("ascii", errors="ignore")).hexdigest()	

def split_var_name(name):
    items = re.match(r"([A-Z]+)([0-9]+)([A-Z]+)", name, re.I).groups()
    return items[0],int(items[1]),items[2]	