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
    
    IsaB=1,
    
    blaZ=2,
    
    lnuA=3,
    
    lnuB=4,
    
    ermB=5,
    
    ermC=6,
    
    ermA=7,
    
    ermY=8,
    
    mupB=9,
    
    mupA=10,
    
    ermT=11,
    
    dfrA=12,
    
    vgaALC=13,
    
    dfrC=14,
    
    dfrD=15,
    
    dfrG=16,
    
    mecC=17,
    
    mecA=18,
    
    dfrK=19,
    
    vgaB=20,
    
    vgaA=21,
    
    tetK=22,
    
    tetM=23,
    
    tetL=24,
    
    tetO=25,
    
    msrA=26,
    
    vanA=27,
    
    vanC=28,
    
    vanB=29,
    
    fusB=30,
    
    fusC=31,
    
    str=32,
    
    aacAaphD=33,
    
    PVL=34,
    
    unspecified_gpg = 35
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 35    //ignore unspecified_gpg
#define MAX_LEN_GENE 3155

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
  char* fasta_id;  
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