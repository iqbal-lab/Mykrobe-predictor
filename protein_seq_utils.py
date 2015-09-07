import itertools
from models import subsitution
from Bio.Seq import Seq
def checkForMismatchSubsitutions(seq,ref):
    subsitution_list = []
    assert len(seq) == len(ref)
    for i,(ref,sample) in enumerate(itertools.izip(ref,seq)):
        if not ref == sample:
            subsitution_list.append(subsitution(ref,sample,i+1)) ## add 1 for base one consistancy    
    print subsitution_list
    return subsitution_list

def sequenceToListTriplets(seq):
    outlist=[]
    tmplist = []
    for i,base in enumerate(list(seq)):
        tmplist.append(base)
        if len(tmplist) == 3 :
            outlist.append(tmplist)
            tmplist = []
    return outlist

def frameShiftCheck(triplet_list):
    """Looks for framshift causing deletions"""
    framshift_pos = []
    for i,triplet in enumerate(triplet_list):
        seq_string = "".join(triplet)
        if "_" in triplet and not seq_string == "___":
            framshift_pos.append(i)
    return framshift_pos


def translateWithINDELS(triplet_list):
    protein_seq = []
    for i,triplet in enumerate(triplet_list):
        assert len(triplet) == 3
        seq_string = "".join(triplet)
        if seq_string == "___":
            protein_seq.append("_")
        elif "_" in triplet and not seq_string == "___":
            logging.info("Framshift at position %i?" % i)
            0/0
        else:
            protein_seq.append(str(Seq(seq_string).translate()))
    return "".join(protein_seq)

