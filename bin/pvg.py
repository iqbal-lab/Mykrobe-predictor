from utils import AutoVivification,kmerCombo,getComboMuts,unique
import Bio
from Bio.Data import CodonTable 
from Bio.SeqRecord import SeqRecord
from Bio import Entrez,SeqIO,Alphabet,Seq
from collections import namedtuple
from copy import copy
import logging
subsitution = namedtuple('Subsitution', 'ref alt pos')
insertion = namedtuple('Insertion', 'ref ins pos')
def make_backward_codon_table():
    table = AutoVivification()
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    for kmer in kmerCombo(3):
        if kmer not in standard_table.stop_codons:
            try:
                table[standard_table.forward_table[kmer]].append(kmer)
            except:
                table[standard_table.forward_table[kmer]] = [kmer]
    return table


def introduce_amino_acid_subsitution(gene,ref,sub,identifier=""):
    """Function to take reference protien sequence, an amino acid subsitution namedtuple"""
    ref_protien_seq =  ref.seq.translate()
    ref_DNA_seq = ref.seq

    prot_list= list(str(ref_protien_seq))
    DNA_list = list(str(ref_DNA_seq))
    prot_list[sub.pos] = sub.alt
    ## Make an equivalent subsitution in the DNA
    backward_codon_table = make_backward_codon_table()

    long_id = "%s_%s" % (gene,sub.ref+str(sub.pos+1)+sub.alt)
    ## Check if the mutation only exists within combo of others
    comboMutDic = getComboMuts()

    if long_id in comboMutDic:
        combo = comboMutDic[long_id]
    else:
        combo = '0'



    codonList = backward_codon_table[sub.alt]
    seqOutList = [SeqRecord(ref_DNA_seq[(sub.pos*3)-30:(sub.pos*3)+33],id="ref_"+identifier+"_sub-%i-%s" % (len(codonList),gene), description="",name = "ref_"+identifier)]
    assert len(seqOutList[0]) == 63
    for j,kmer in enumerate(codonList):
    ## Go through the possible codons for the alt aminoacid and generate sequences
        tmp_DNA_list = DNA_list
        tmp_DNA_list[sub.pos*3 : sub.pos*3 +3] = list(kmer)

        seqOut = SeqRecord(Seq.Seq("".join ( tmp_DNA_list),Alphabet.IUPAC.unambiguous_dna),id= "".join(kmer)+'_'+identifier+"_alt-%i-%s" %((j+1),gene ) ,description="",name = "".join(kmer)+'_'+identifier)
        # try:
        seqOutList.append(seqOut[(sub.pos*3)-30:(sub.pos*3)+33])
        assert len(seqOutList[-1]) == 63 
        # except:
        #     diff =  len(tmp_DNA_list)-sub.pos*3
        #     newSeq = seqOut[(sub.pos*3)-30+diff-33:(sub.pos*3)+33]
        #     if  len(newSeq) == 63:
        #         seqOutList.append(newSeq)
        ## Make sure the new sequence tranlates to the alternate 
        
        assert str(seqOut.seq.translate()) == "".join(prot_list)
        assert str(seqOut.seq.translate()) != str(ref_protien_seq)

    return seqOutList

def introduce_amino_acid_insertion(gene,ref,sub,identifier=""):
    """Function to take reference protien sequence, an amino acid subsitution namedtuple"""
    ref_protien_seq =  ref.seq.translate()
    ref_DNA_seq = ref.seq

    prot_list= list(str(ref_protien_seq))
    DNA_list = list(str(ref_DNA_seq))
    prot_list.insert(sub.pos,sub.ins)
    ## Make an equivalent subsitution in the DNA
    backward_codon_table = make_backward_codon_table()

    long_id = "%s_%s" % (gene,"ins"+str(sub.pos+1)+sub.ins)
    ## Check if the mutation only exists within combo of others
    comboMutDic = getComboMuts()

    if long_id in comboMutDic:
        combo = comboMutDic[long_id]
    else:
        combo = '0'

    codonList = backward_codon_table[sub.ins]
    seqOutList = [SeqRecord(ref_DNA_seq[(sub.pos*3)-30:(sub.pos*3)+33],id="ref_"+identifier+"_ins-%i-%s" % (len(codonList),gene), description="",name = "ref_"+identifier)]
    assert len(seqOutList[0]) == 63
    for j,kmer in enumerate(codonList):
    ## Go through the possible codons for the alt aminoacid and generate sequences
        tmp_DNA_list = list(str(ref_DNA_seq))
        for i,base in enumerate(kmer):
            tmp_DNA_list.insert(sub.pos*3+i,base)

        seqOut = SeqRecord(Seq.Seq("".join ( tmp_DNA_list),Alphabet.IUPAC.unambiguous_dna),id = "".join(kmer)+'_'+identifier+"_alt-%i-%s" % ((j+1),gene),description="",name = "".join(kmer)+'_'+identifier)
        seqOutList.append(seqOut[(sub.pos*3)-30:(sub.pos*3)+36])
        ## Make sure the new sequence tranlates to the alternate 
        assert len(seqOutList[-1]) == 66
        assert str(seqOut.seq.translate()) == "".join(prot_list)
        assert str(seqOut.seq.translate()) != str(ref_protien_seq)

    return seqOutList


def introduce_amino_acid_change(gene,ref,sub,identifier=""):
    if type(sub).__name__ == "Subsitution":
        return introduce_amino_acid_subsitution(gene,ref,sub,identifier)
    elif type(sub).__name__ == "Insertion":
        return introduce_amino_acid_insertion(gene,ref,sub,identifier)
    else:
        0/0





    

class Panel(object):
    """A panel consists of multiple drug panels"""
    def __init__(self,organism,drugPanelList):
        self.organism = organism
        self.drugPanelList = drugPanelList

    @property 
    def geneList(self):
        geneList = []
        for drugPanel in self.drugPanelList:
            for variant in drugPanel.variantList:
                geneList.append(variant.gene.name)
        return unique(geneList)


    def writeAllPanelFastas(self,outdir):
        for drugPanel in self.drugPanelList:
          SeqIO.write(drugPanel.getSeqRecordList(),outdir+"%s.fa"%drugPanel.drug,'fasta')

    def writeDrugEmum(self,outdir):
        drugEnums = []
        for i,drugPanel in enumerate(self.drugPanelList):
            if not i == len(self.drugPanelList)-1:
                drugEnum = "%s=%i,"% (drugPanel.drug,i+1)
            else:
                drugEnum = "%s=%i"% (drugPanel.drug,i+1)
            drugEnums.append(drugEnum )

        base = """typedef enum 
                 {
                   NoDrug=0,
                   %s
                  } Antibiotic;
                  """ % "\n".join(drugEnums)
        with open(outdir+'antibioticEnum.c','w') as outfile:
            outfile.write(base)

    def write_map_antibiotic_enum_to_str(self,outdir):
        elseifList = []
        for i,drugPanel in enumerate(self.drugPanelList):
            elseif = """else if (ab==%s)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "%s");
                }
        """ % (drugPanel.drug,drugPanel.drug)
            elseifList.append(elseif)
        


        base = """
        void map_antibiotic_enum_to_str(Antibiotic ab, StrBuf* name)
{
        if (ab==NoDrug)
        {
          strbuf_reset(name);
          strbuf_append_str(name, "NoDrug");
        }
        %s 
        else
        {
          die("Impossible - compiler should not allow this\n");
        }
}

        """ % "\n\t".join(elseifList)
        with open(outdir+"map_antibiotic_enum_to_str.c",'w') as outfile:
            outfile.write(base)


    def write_is_drug_susceptible(self,outdir):
        for drugPanel in self.drugPanelList:
            isDrugSusceptible = """
    Troolean is_%s_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = %s;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/%s/antibiotics/%s.fa");

      abi->num_mutations = %i;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = %s;
      int last_mut = %s;

      int i;
       //if you have any of these resistance alleles - call resistant
  double max_conf=-9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (max_conf < best_model.conf)
    {
      max_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (max_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}""" % (drugPanel.drug, drugPanel.drug, self.organism, drugPanel.drug, 
                drugPanel.num_mutations, drugPanel.first_mut, 
                drugPanel.last_mut)

            isDrugSusceptibleh = """Troolean is_%s_susceptible(dBGraph* db_graph,
                     int (*file_reader)(FILE * fp, 
                            Sequence * seq, 
                            int max_read_length, 
                            boolean new_entry, 
                            boolean * full_entry),
                     ReadingUtils* rutils,
                     ResVarInfo* tmp_rvi,
                     GeneInfo* tmp_gi,
                     AntibioticInfo* abi,
                     StrBuf* install_dir,
                     int ignore_first, int ignore_last, int expected_covg,
                     double lambda_g, double lambda_e, double err_rate
                     );"""  % drugPanel.drug
            with open(outdir + "is_%s_susceptible.c" % drugPanel.drug , 'w' ) as outfile:
                outfile.write(isDrugSusceptible)
            with open(outdir + "is_%s_susceptible.h" % drugPanel.drug , 'w' ) as outfile:
                outfile.write(isDrugSusceptibleh)
    def write_GeneMutationGene_Enum(self,outdir):
        geneEnumList = []
        geneifElseList = []
        for i,gene in enumerate(self.geneList):
            if i == len(self.geneList)-1 :
                geneEnumList.append("%s = %i" % (gene,i+1) )
            else:
                geneEnumList.append("%s = %i," % (gene,i+1) )

            if i == 0 :
                geneifElse = """if (strcmp(name->buff, "%s")==0)
                    {
                      return %s ;
                    }""" % (gene,gene)
            else:
                geneifElse = """
                else if (strcmp(name->buff, "%s")==0)
                {
                  return %s;
                }
                """ % (gene,gene)
            geneifElseList.append(geneifElse)



        baseEnum = """typedef enum
          {
            Unknown = 0,
            %s
            
          }GeneMutationGene;""" % "\n\t".join(geneEnumList)

        ifelse = """GeneMutationGene map_gene_name_str_to_genename(StrBuf* name)
            {
            %s
                else
                {
                  die("Failed to parse gene name - got %s", name->buff);
                }

            }
            """ % ("\n\t".join(geneifElseList),"%s")
        with open(outdir+"GeneMutationGene_Enum.h",'w') as outfile:
            outfile.write(baseEnum)
        with open(outdir + "map_gene_name_str_to_genename.c",'w') as outfile:
            outfile.write(ifelse)

    def write_KnownMutation_enum(self,outdir):
        enumList = []
        mutifElseList = []
        variantidentifierList = []
        cnt = 0
        for drugPanel in self.drugPanelList:
            for variant in drugPanel.variantList:
                if not variant.identifier in variantidentifierList:
                    variantidentifierList.append(variant.identifier)
                    enum = "%s_%s = %i," %(variant.gene.name, variant.identifier, cnt)
                    cnt+=1
                    enumList.append(enum)
                    if cnt == 1 :
                        mutifElse = """
                        if ( (strcmp(sbuf->buff, "%s")==0) && (gene==%s) )
                          {
                            return %s_%s;
                          }
                        """ % (variant.identifier,variant.gene.name,
                                variant.gene.name,variant.identifier)
                    else:
                        mutifElse = """
                        else if ( (strcmp(sbuf->buff, "%s")==0) && (gene==%s) )
                          {
                            return %s_%s;
                          }
                        """ % (variant.identifier,variant.gene.name,
                                variant.gene.name,variant.identifier)
                    mutifElseList.append(mutifElse)

        baseKnownMutEnum = """typedef enum
    {
    %s
    NotSpecified = %i,
  } KnownMutation;
  #define NUM_KNOWN_MUTATIONS %i """ % ("\n\t".join(enumList),len(enumList),len(enumList))

        ifelse = """KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene)
{
%s
                    else 
                        {
                          die("Parsing error - unknown mutation %s", sbuf->buff);
                          return NotSpecified;
                        } 

}""" % ("\n\t".join(mutifElseList),"%s")
        with open(outdir+'KnownMutation_enum.h','w') as outfile:
            outfile.write(baseKnownMutEnum)
        with open(outdir+"map_mutation_name_to_enum.c",'w') as outfile:
            outfile.write(ifelse)

    def write_main_print_drug_susc(self,outdir):
        mainList = []
        for drugPanel in self.drugPanelList:
            base = """print_antibiotic_susceptibility(db_graph, &file_reader_fasta, ru, tmp_rvi, tmp_gi, abi,
                      &is_%s_susceptible, tmp_name, cmd_line->install_dir,
                      ignore, ignore, expected_depth, lambda_g_err, lambda_e_err, err_rate, cmd_line->format, output_last); """ % (drugPanel.drug)
            mainList.append(base)
        with open(outdir + "main.c",'w') as outfile:
            outfile.write("\n".join(mainList))








class DrugPanel(object):
    """A drug panel consists of several resistance variants"""
    def __init__(self,drug,resistanceVariantList):
        self.drug = drug
        self.variantList = resistanceVariantList

    @property 
    def orderedVariantList(self):
      """Sort the variants so that all variants with the same identifier are together"""
      orderedVariantList = []
      for identifier in self.variant_identifiers:
        tmpVariantList = [v for v in self.variantList if v.identifier == identifier]
        ## Reduce variant list 
        tmpRefList = [str(v.getReferenceAllele().seq) for v in tmpVariantList]
        ## Get the indexes of the unique references 
        uniqueRefs = unique(tmpRefList)
        uniqueRefIndexes = [tmpRefList.index(u) for u in uniqueRefs]
        ## Use these indexes to filter the variants
        tmpVariantList = [tmpVariantList[i] for i in uniqueRefIndexes]
        orderedVariantList.extend(tmpVariantList)
      return orderedVariantList


    def getSeqRecordList(self):
        logging.debug("%s"%self.drug)
        drugfastalist = []
        ### The sequences must be ordered such that all of the variants with the same IDs are together. 
        for variant in self.orderedVariantList:
            logging.debug("%s" % " ".join([variant.identifier, variant.gene.name,variant.gene.ref_seq_id]))
            try:
              referenceAllele = variant.getReferenceAllele()
              alternateAlleles = variant.getAlternateAlleles()
              if (referenceAllele is not None) and (alternateAlleles is not None):
                  drugfastalist.append(referenceAllele)
                  drugfastalist.extend(alternateAlleles)
              else:
                ### If there's no reference or alt alleles it is because the reference is one of the resistant amino acids. 
                ### This means that the "referenceAllele" is resistance. 
                ### We need to generate all the possible back mutations which would give the S allele.
                logging.debug("%s" % " ".join(["else",variant.identifier, variant.gene.name,variant.gene.ref_seq_id]))
                referenceAlleles = variant.getReferenceAllelesIfRefIsResistance()
                alternateAllele = variant.getAlternateAlleleIfRefIsResistance()
                for referenceAllele in referenceAlleles:
                    drugfastalist.append(referenceAllele)
                    drugfastalist.append(alternateAllele)
            except:
              pass

        return drugfastalist

    @property 
    def variant_identifiers(self):
      return unique([v.identifier for v in self.variantList])

    @property
    def num_mutations(self):
        return len(self.variant_identifiers)

    @property 
    def first_mut(self):
        firstVariant = self.variantList[0]
        return "%s_%s" % (firstVariant.gene.name,firstVariant.identifier )

    @property
    def last_mut(self):
        lastVariant = self.variantList[-1]
        return "%s_%s" % (lastVariant.gene.name,lastVariant.identifier )



class ComboResistanceVariants(object):
    """N variants which together cause induces resistance"""
    def __init__(self,resistanceVariantList,kmer):
        self.kmer = kmer
        self.rvl = resistanceVariantList
        self.rv1 = self.rvl [0]
        self.rv2 = self.rvl [1]
        assert self.rv1.gene == self.rv2.gene
        self.gene = self.rv1.gene
        assert self.rv1.pos < self.rv2.pos
        self.alphabet = resistanceVariantList[0].alphabet
        self.identifier = self.rv1.identifier+"and"+self.rv2.identifier
        ## Asserts
        assert len(self.rvl) == 2
        assert self.rv1 != self.rv2
        if self.alphabet == Alphabet.IUPAC.protein:
            ## The variants must be overlapping in dna space
            assert abs(self.rv1.pos - self.rv2.pos)+1 <= (self.kmer/3) 
        elif self.alphabet == Alphabet.IUPAC.unambiguous_dna:
            assert abs(self.rv1.pos - self.rv2.pos)+1 <= self.kmer


    def getReferenceAllele(self,l=63):
        s =  self.rvl[0].getReferenceAllele(l=l)
        s2 = self.rvl[0].getReferenceAllele(l=l)
        if s is not None:
          s.id="ref_"+self.identifier+"_ref-%i-%s" % ( (self.rv1.num_alternates*self.rv2.num_alternates),self.gene.name)
        if s is None or s2 is None:
          return None
        else:
          return s

    def getReferenceAllelesIfRefIsResistance(self,l=63):
        # return self.rvl[0].getReferenceAllelesIfRefIsResistance(l=l)
        print ("Skipping combo mutation where one of the resistant alleles is in the reference")
        return Seq.Seq("")

    def getAlternateAlleleIfRefIsResistance(self,l=63):
        # return self.rvl[0].getAlternateAlleleIfRefIsResistance(l=l)
        print ("Skipping combo mutation where one of the resistant alleles is in the reference")
        return Seq.Seq("")

    def getAlternateAlleles(self,l=63):
        """Induce all of the variants simultaniously"""
        if self.alphabet == Alphabet.IUPAC.unambiguous_dna:
            reftmp = self.getReferenceAllele()
            if reftmp is not None:
              ref = list(reftmp.seq)
              alternateAlleles = []
              cnt=0
              for alt1 in self.rv1.alternates:
                  for alt2 in self.rv2.alternates:
                      cnt+=1
                      alternateAllele = copy(ref)
                      alternateAllele[ ((l+1)/2)-1] = alt1
                      alternateAllele[((l+1)/2)-1 + (self.rv2.pos - self.rv1.pos)] = alt2
                      assert alternateAllele != ref
                      alternateAlleles.append(  SeqRecord( Seq.Seq("".join(alternateAllele)),
                                                   id="alt_"+self.identifier+"_alt-%i-%s" % (cnt,self.rv1.gene.name), name="",
                                                      description="")
                                                )
              return alternateAlleles
            else:
              return None

        elif self.alphabet == Alphabet.IUPAC.protein:
            ### 32+- midel position of the codon
            alternateAlleles = []
            ## For each reference
            ref = self.getReferenceAllele(l=l)
            if ref is not None:
              ref = list(self.getReferenceAllele().seq)
              cnt=0
              for alt1 in self.rv1.alternates:
                  codonList1 = self.rv1.backward_codon_table[alt1]
                  for alt2 in self.rv2.alternates:
                      ## get the codons for this alternate
                      codonList2 = self.rv2.backward_codon_table[alt2]
                      for codon1 in codonList1:
                          for codon2 in codonList2:
                              alternateAllele = copy(ref)
                              ## Replace the 3 middle based with the new codon
                              listRange1 = (  ((l+1)/2) -2, (l+1)/2+1 )
                              listRange2 = listRange1[0] + (self.rv2.pos - self.rv1.pos)*3, listRange1[1]+(self.rv2.pos - self.rv1.pos)*3

                              replaceTrimer1 = "".join(ref[listRange1[0]:listRange1[1]])
                              replaceTrimer2 = "".join(ref[listRange2[0]:listRange2[1]])

                              alternateAllele[listRange1[0]:listRange1[1]] = codon1
                              alternateAllele[listRange2[0]:listRange2[1]] = codon2

                              assert alternateAllele != ref
                              alternateAllele = Seq.Seq("".join(alternateAllele))
                              assert alternateAllele not in alternateAlleles
                              cnt+=1
                              alternateAlleles.append(SeqRecord ( alternateAllele,
                                                              id="alt_"+self.identifier+"_alt-%i-%s" % (cnt,self.rv1.gene.name), name="",
                                                                  description="" ) )
            else:
                alternateAlleles = None
            return alternateAlleles


        


class Gene(object):
    """A gene in which amino acid mutations occur"""
    def __init__(self,name,ref_seq,ref_seq_id,start_pos,end_pos,strand="forward"):
        self.name = name
        self.ref_seq = ref_seq
        self.ref_seq_id = ref_seq_id
        ## Start position and end position are given in 1-based coords.
        self.start_pos = int(start_pos)-1
        self.end_pos = int(end_pos)
        self.strand = strand
        if strand == "forward":
            self.forward = True
        elif strand == "reverse":
            self.forward = False
        else:
            raise ValueError("strand must be either forward or reverse")

    @property 
    def gene_length(self):
      return self.end_pos-self.start_pos

    @property 
    def dna(self):
        """The dna sequence of the gene"""
        if self.forward:
            return self.ref_seq[self.start_pos:self.end_pos]
        else:
            return (self.ref_seq[self.start_pos:self.end_pos]).reverse_complement()
    @property 
    def prot(self):
        """The protien sequence of the gene"""
        try:
            if "N" in self.dna[-3:] or "N" in self.dna[:3]:
                return self.dna.translate(cds=False,table=11)
            else:
                return self.dna.translate(cds=True,table=11)
        except Bio.Data.CodonTable.TranslationError,e :
            if self.forward:
                # print self.ref_seq[self.start_pos:self.end_pos]
                raise e
            else:
                # print (self.ref_seq[self.start_pos:self.end_pos]).reverse_complement()
                raise e
    def has_position(self,pos):
        return pos <= self.end_pos and pos >= self.start_pos



    def getBase(self,pos):
        """returns the DNA base at a particular position given in 1 based"""
        if pos > 0:
            if pos <= self.gene_length:
                return self.dna[pos-1]
            else:
                ## If we're outside the length of the gene we go back to the reference
                return self.ref_seq[self.start_pos+pos]
        elif pos < 0:
            """If we get a negative position we need to go back to the reference"""
            if self.forward:
                return self.ref_seq[self.start_pos+pos]
            else:
                ## If we're looking at the reverse strand then upstream is relative to the endposition. 
                ## I.e gene -10 is 10 bases further then end position.
                ## This means that if end position is say, 5332, then we want position 5342 in 1-based. 5341 in 0-based
                return self.ref_seq[self.end_pos - pos -1]
        elif pos == 0:
            """Not sure if 0 positions in this context make sense: pos 1 is the start of the gene. -1 is 1 upstream."""
            raise(IndexError("Positions are 1 based - pos = 0 doesn't make sense here. "))
    def getAminoAcid(self,pos):
        """Returns the aminoacid at a particular 1 based potision"""
        return self.prot[pos-1]

    def getDNARange(self,pos,l):
        positions  = range(pos-((l-1)/2),pos+((l-1)/2)+1)
        if 0 in positions:
            positions.append(positions[-1]+1)
            positions.remove(0)

        string = "".join([self.getBase(p) for p in positions])
        assert len(string) == l
        return Seq.Seq(string)

class ResistanceVariant(object):
    """A variant in a reference sequence which induces resistance to a drug"""
    def __init__(self,drug,gene,alphabet,pos,ref,alts):
        self.ref_seq = gene.ref_seq
        self.drug = drug.strip()
        self.gene = gene
        self.pos = int(pos)
        self.ref = ref
        self.alts = alts
        self.identifier = "%s%i%s" % (self.ref,self.pos,"".join(self.alts))
        self.identifier = self.identifier.replace("*","X")
        self.identifier = self.identifier.replace("-","u")

        assert ref not in alts
        if not self.ref and self.alts:
          raise NotImplementedError("No rule for INDELs yet.")
        if alphabet =="DNA" or alphabet == Alphabet.IUPAC.unambiguous_dna:
            self.alphabet = Alphabet.IUPAC.unambiguous_dna
            self.gene_reference = gene.getBase(self.pos)
            if self.alts == "*":
                self.alternates = ['A','T','C','G']
                self.alternates.remove(self.ref)
                assert len(self.alternates) == 3
            else:
                self.alternates = self.alts
                
        elif alphabet == "PROTEIN" or alphabet == Alphabet.IUPAC.protein:
            self.alphabet = Alphabet.IUPAC.protein
            self.gene_reference = gene.getAminoAcid(self.pos)
            self.backward_codon_table = make_backward_codon_table()
            if self.alts == "*":
                self.alternates = ['A','R','N','D','C','E','Q','G'\
                                    ,'H','I','L','K','M','F','P','S','T',\
                                    'W','Y','V']
                self.alternates.remove(self.ref)
                assert len(self.alternates) == 19
            else:
                self.alternates = self.alts
        else:
            raise ValueError("alphabet must be either DNA or PROTEIN")
        ## Ensure that the positions make sense
        # if not (gene_reference == self.ref or gene_reference in ["X"]):

        # 
    ## Get the dna position
    def getDNAPos(self):
        return (self.gene.start_pos + (self.pos-1)*3,self.gene.start_pos + (self.pos)*3)
    def getCodon(self):
        codonRange = self.getDNAPos()
        return self.gene.ref_seq[codonRange[0]:codonRange[1]]

    def getReferenceAllele(self,l=63):
        ## Assert that the reference allele is not that alternate
        if self.gene_reference in self.alts:
          logging.info("Reference codon %s is one of the  resistant aminoacids %s" % (self.gene_reference,self.alts))
          ## When the reference codon is the resitance allele generate the back mutation to create the S alelle. 
          # rv = ResistanceVariant(drug=self.drug,gene=self.gene,alphabet=self.alphabet,pos=self.pos,ref=self.alt,alt=self.ref)
          ## Generate all of the S alleles, the alternate alleles in the back mutation
          return None
        else:
          # assert (self.gene_reference == self.ref or self.gene_reference in ["X"])
          if self.alphabet == Alphabet.IUPAC.unambiguous_dna:
              s = SeqRecord ( self.gene.getDNARange(pos=self.pos,l=l) , 
                                  id="ref_"+self.identifier+"_sub-%i-%s" % (self.num_alternates,self.gene.name), name="",
                                  description="")
              assert not "N" in str(s.seq) ## No N's allowed
              return s
          elif self.alphabet == Alphabet.IUPAC.protein:
              ### 32+- midel position of the codon
              pos = self.pos * 3 -1
              s =  SeqRecord ( self.gene.getDNARange(pos=pos,l=l) ,
                                  id="ref_"+self.identifier+"_sub-%i-%s" % (self.num_alternates,self.gene.name), name="",
                                  description="")

              assert not "N" in str(s.seq) ## No N's allowed
              return s


    def getAlternateAlleleIfRefIsResistance(self,l=63):
        ### Return the alelle of the reference genome ensuring that it carries the resistant allele and that the 
        ### fasta description is of the alternate of a subsitution with one alternate allele.
        assert self.gene_reference in self.alts
        if self.alphabet == Alphabet.IUPAC.unambiguous_dna:
          s = SeqRecord ( self.gene.getDNARange(pos=self.pos,l=l) , 
                                  id="alt_"+self.identifier+"_alt-%i-%s" % (1,self.gene.name), name="",
                                  description="")
        elif self.alphabet == Alphabet.IUPAC.protein:
          pos = self.pos * 3 -1
          s =  SeqRecord ( self.gene.getDNARange(pos=pos,l=l) ,
                                  id="alt_"+self.identifier+"_alt-%i-%s" % (1,self.gene.name), name="",
                                  description="")
        assert not "N" in str(s.seq) ## No N's allowed
        return s
    def getReferenceAllelesIfRefIsResistance(self,l=63):
        ### Return the putative susceptible alleles given that the reference carries the resistant
        assert self.gene_reference in self.alts
        if self.alphabet == Alphabet.IUPAC.unambiguous_dna:
                logging.warning("Haven't written DNA version of this case")
                return 0/0
        elif self.alphabet == Alphabet.IUPAC.protein:
            ### 32+- midel position of the codon
            referenceAlleles = []
            ref = list(self.getAlternateAlleleIfRefIsResistance().seq)
            alt = self.ref ## The susceptible allele
            ## get the codons for this alternate
            codonList = self.backward_codon_table[alt]
            for codon in codonList:
                alternateAllele = copy(ref)
                ## Replace the 3 middle based with the new codon
                listRange = (  ((l+1)/2) -2, (l+1)/2+1 )
                replaceTrimer = "".join(ref[listRange[0]:listRange[1]])
                # assert replaceTrimer in self.backward_codon_table[self.ref]
                if not replaceTrimer in self.backward_codon_table[self.ref]:
                    logging.warning("Reference trimer not in backward codon list")
                    logging.warning("%s" % [replaceTrimer,self.backward_codon_table[self.ref]])
                alternateAllele[listRange[0]:listRange[1]] = codon
                assert alternateAllele != ref
                alternateAllele = Seq.Seq("".join(alternateAllele))
                assert alternateAllele not in referenceAlleles
                ## Each of the references have one alternate (the resistant allele)
                referenceAlleles.append(SeqRecord ( alternateAllele,
                                                id="ref_"+self.identifier+"_ref-%i-%s" % (1,self.gene.name), name="",
                                                    description="" ) )

        return referenceAlleles

            

    def getAlternateAlleles(self,l=63):
        if self.gene_reference in self.alts:
            print "Reference codon %s is one of the resistant aminoacid %s" % (self.gene_reference,self.alts)
            return None
        else:
            if self.alphabet == Alphabet.IUPAC.unambiguous_dna:
                ref = list(self.getReferenceAllele().seq)
                alternateAlleles = []
                cnt=0
                for alt in self.alternates:
                    cnt+=1
                    alternateAllele = copy(ref)
                    alternateAllele[ ((l+1)/2)-1] = alt
                    assert alternateAllele != ref
                    alternateAlleles.append(  SeqRecord( Seq.Seq("".join(alternateAllele)),
                                                 id="alt_"+self.identifier+"_alt-%i-%s" % (cnt,self.gene.name), name="",
                                                    description="")
                                              )
                return alternateAlleles
            elif self.alphabet == Alphabet.IUPAC.protein:
                ### 32+- midel position of the codon
                alternateAlleles = []
                ref = list(self.getReferenceAllele().seq)
                cnt=0
                for alt in self.alternates:

                    ## get the codons for this alternate
                    codonList = self.backward_codon_table[alt]
                    for codon in codonList:
                        alternateAllele = copy(ref)
                        ## Replace the 3 middle based with the new codon
                        listRange = (  ((l+1)/2) -2, (l+1)/2+1 )
                        replaceTrimer = "".join(ref[listRange[0]:listRange[1]])
                        # assert replaceTrimer in self.backward_codon_table[self.ref]
                        # if not replaceTrimer in self.backward_codon_table[self.ref]:
                        #     print"Reference trimer not in backward codon list" 
                        #     print replaceTrimer,self.backward_codon_table[self.ref]
                        alternateAllele[listRange[0]:listRange[1]] = codon
                        assert alternateAllele != ref
                        alternateAllele = Seq.Seq("".join(alternateAllele))
                        assert alternateAllele not in alternateAlleles
                        cnt+=1
                        alternateAlleles.append(SeqRecord ( alternateAllele,
                                                        id="alt_"+self.identifier+"_alt-%i-%s" % (cnt,self.gene.name), name="",
                                                            description="" ) )
                return alternateAlleles

    @property 
    def num_alternates(self):
        """Calculates the number of alternate alelles"""
        if self.alphabet == Alphabet.IUPAC.unambiguous_dna:
            return len(self.alternates)
        elif self.alphabet == Alphabet.IUPAC.protein:
            totalCodons = []
            for alt in self.alternates:
                totalCodons.extend(self.backward_codon_table[alt])
            return len(totalCodons)



        
        



