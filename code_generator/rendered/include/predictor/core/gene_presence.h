/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  gene_presence.h - build the de Bruijn graph
*/

#ifndef GENE_PRESENCE_H
#define GENE_PRESENCE_H


#include "global.h" // Covg def 
#include "seq.h" // Sequence def 
#include "binary_kmer.h" // KmerSlidingWindow
#include "dB_graph.h"// dBGraph
#include "dB_graph_supernode.h"




typedef enum
  {
    
    vgbA=0,
    
    aacAaphD=1,
    
    IsaB=2,
    
    aadEant6Ia=3,
    
    aadDaph4Ia=4,
    
    qacB=5,
    
    qacA=6,
    
    cfr=7,
    
    blaZ=8,
    
    lnuA=9,
    
    lnuB=10,
    
    ermB=11,
    
    ermC=12,
    
    ermA=13,
    
    ermY=14,
    
    mupB=15,
    
    mupA=16,
    
    aph2Ic=17,
    
    ant9Ia=18,
    
    mecA=19,
    
    ermT=20,
    
    ant9Ib=21,
    
    aphA3aph3III=22,
    
    dfrA=23,
    
    vgaALC=24,
    
    dfrC=25,
    
    dfrD=26,
    
    dfrG=27,
    
    mecC=28,
    
    aad9spc=29,
    
    dfrK=30,
    
    vgaB=31,
    
    vgaA=32,
    
    tetK=33,
    
    tetM=34,
    
    tetL=35,
    
    tetO=36,
    
    msrA=37,
    
    mphC=38,
    
    sat4=39,
    
    cat=40,
    
    vanA=41,
    
    vanC=42,
    
    vanB=43,
    
    fusB=44,
    
    fusC=45,
    
    str=46,
    
    qacCsmr=47,
    
    luk=48,
    
    unspecified_gpg = 49
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 49    //ignore unspecified_gpg
#define MAX_LEN_GENE 3110

GenePresenceGene map_string_to_gene_presence_gene(StrBuf* sbuf);
boolean map_gene_to_fasta(GenePresenceGene gene, StrBuf* fa, StrBuf* install_dir);

typedef struct
{
  Covg median_covg;
  Covg min_covg;
  Covg median_covg_on_nonzero_nodes;
  int num_gaps;
  int len;
  int  percent_nonzero;
  StrBuf* strbuf;
  GenePresenceGene name;
} GeneInfo;

GeneInfo* alloc_and_init_gene_info();
void free_gene_info(GeneInfo* rvi);
void reset_gene_info(GeneInfo* rvi);
void copy_gene_info(GeneInfo* from_gi, GeneInfo* to_gi);

int get_next_gene_info(FILE* fp, 
		       dBGraph* db_graph, 
		       GeneInfo* gene_info,
		       Sequence* seq, 
		       KmerSlidingWindow* kmer_window,
		       int (*file_reader)(FILE * fp, 
					  Sequence * seq, 
					  int max_read_length, 
					  boolean new_entry, 
					  boolean * full_entry),
		       dBNode** array_nodes, 
		       Orientation*  array_or,
		       CovgArray* working_ca, int max_read_length);
const char* map_enum_to_gene_name(GenePresenceGene gene);


#endif