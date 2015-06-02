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
    
    qacB=3,
    
    qacA=4,
    
    cfr=5,
    
    blaZ=6,
    
    lnuA=7,
    
    lnuB=8,
    
    ermB=9,
    
    ermC=10,
    
    ermA=11,
    
    ermY=12,
    
    mupB=13,
    
    mupA=14,
    
    ermT=15,
    
    ant9Ib=16,
    
    dfrA=17,
    
    vgaALC=18,
    
    dfrC=19,
    
    dfrD=20,
    
    dfrG=21,
    
    mecC=22,
    
    mecA=23,
    
    dfrK=24,
    
    vgaB=25,
    
    vgaA=26,
    
    tetK=27,
    
    tetM=28,
    
    tetL=29,
    
    tetO=30,
    
    msrA=31,
    
    sat4=32,
    
    cat=33,
    
    vanA=34,
    
    vanC=35,
    
    vanB=36,
    
    fusB=37,
    
    fusC=38,
    
    str=39,
    
    qacCsmr=40,
    
    arcA=41,
    
    arcB=42,
    
    arcC=43,
    
    arcD=44,
    
    ccrA=45,
    
    ccrB=46,
    
    ccrCa=47,
    
    ccrCb=48,
    
    ccrCc=49,
    
    chp=50,
    
    eta=51,
    
    etb=52,
    
    etd=53,
    
    luk=54,
    
    lukM=55,
    
    lukMF=56,
    
    lukPVF=57,
    
    lukPVS=58,
    
    sak=59,
    
    sasX=60,
    
    scn=61,
    
    sea=62,
    
    seb=63,
    
    sec=64,
    
    sed=65,
    
    see=66,
    
    seg=67,
    
    seh=68,
    
    sei=69,
    
    sej=70,
    
    selR=71,
    
    sep=72,
    
    seu=73,
    
    tsst1=74,
    
    unspecified_gpg = 75
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 75    //ignore unspecified_gpg
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