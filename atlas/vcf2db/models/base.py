import hashlib
def make_hash(s):
	return hashlib.sha256(s.encode("ascii", errors="ignore")).hexdigest()	