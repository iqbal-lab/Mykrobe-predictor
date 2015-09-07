"""DNA variant generator"""
from collections import namedtuple
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet,Seq
insertion = namedtuple('Insertion', 'ins pos')
deletion = namedtuple('Deletion', 'ref pos')


def introduce_dna_insertion(gene,ref,ins,identifier=""):
    """Function to take reference protien sequence, an amino acid subsitution namedtuple"""
    ref_DNA_seq = ref.seq
    DNA_list = list(str(ref_DNA_seq))
    ## Make an equivalent subsitution in the DNA
    long_id = "%s_%s" % (gene,"ins"+str(ins.pos+1)+ins.ins)
    seqOutList = [SeqRecord(ref_DNA_seq[(ins.pos)-30:(ins.pos)+30],id="ref_"+identifier+"_ins--%i--%s--%s" % (1,gene,ins.ins), name="ref_"+identifier)]
    tmp_DNA_list = list(str(ref_DNA_seq))
    tmp_DNA_list.insert(ins.pos,ins.ins)
    seqOut = SeqRecord(Seq.Seq("".join ( tmp_DNA_list),Alphabet.IUPAC.unambiguous_dna),id = "".join(ins.ins)+'_'+identifier+"_alt--%i--%s--%s" % ((0+1),gene,ins.ins),name = "".join(ins.ins)+'_'+identifier)
    seqOutList.append(seqOut[(ins.pos)-30:(ins.pos)+30 + len(ins.ins)])
    assert len(seqOutList[-1].seq) == len(seqOutList[0].seq) + len(ins.ins)
    return seqOutList

def introduce_dna_deletion(gene,ref,deletion,identifier=""):
    """Function to take reference protien sequence, an amino acid subsitution namedtuple"""
    ref_DNA_seq = ref.seq
    DNA_list = list(str(ref_DNA_seq))
    ## Make an equivalent subsitution in the DNA
    long_id = "%s_%s" % (gene,"del"+deletion.ref+str(deletion.pos+1))
    seqOutList = [SeqRecord(ref_DNA_seq[(deletion.pos)-30:(deletion.pos)+30],id="ref_"+identifier+"_del--%i--%s--%s" % (1,gene,deletion.ref), name="ref_"+identifier)]
    tmp_DNA_list = list(str(ref_DNA_seq))
    del tmp_DNA_list[deletion.pos]
    seqOut = SeqRecord(Seq.Seq("".join ( tmp_DNA_list),Alphabet.IUPAC.unambiguous_dna),
    				id = "".join(deletion.ref)+'_'+identifier+"_alt--%i--%s--%s" % ((0+1),
    					gene,deletion.ref),name = "".join(deletion.ref)+'_'+identifier)
    seqOutList.append(seqOut[(deletion.pos)-30:(deletion.pos)+30 - len(deletion.ref)])
    assert len(seqOutList[-1].seq) == len(seqOutList[0].seq) - len(deletion.ref)
    return seqOutList
