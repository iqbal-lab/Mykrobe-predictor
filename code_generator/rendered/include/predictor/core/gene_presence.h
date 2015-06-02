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
    
    aadDaph4Ia=3,
    
    qacB=4,
    
    qacA=5,
    
    cfr=6,
    
    blaZ=7,
    
    lnuA=8,
    
    lnuB=9,
    
    ermB=10,
    
    ermC=11,
    
    ermA=12,
    
    ermY=13,
    
    mupB=14,
    
    mupA=15,
    
    ermT=16,
    
    ant9Ib=17,
    
    dfrA=18,
    
    vgaALC=19,
    
    dfrC=20,
    
    dfrD=21,
    
    dfrG=22,
    
    mecC=23,
    
    mecA=24,
    
    dfrK=25,
    
    vgaB=26,
    
    vgaA=27,
    
    tetK=28,
    
    tetM=29,
    
    tetL=30,
    
    tetO=31,
    
    msrA=32,
    
    mphC=33,
    
    sat4=34,
    
    cat=35,
    
    vanA=36,
    
    vanC=37,
    
    vanB=38,
    
    fusB=39,
    
    fusC=40,
    
    str=41,
    
    qacCsmr=42,
    
    arcA=43,
    
    arcB=44,
    
    arcC=45,
    
    arcD=46,
    
    ccrA=47,
    
    ccrB=48,
    
    ccrCa=49,
    
    ccrCb=50,
    
    ccrCc=51,
    
    chp=52,
    
    eta=53,
    
    etb=54,
    
    etd=55,
    
    luk=56,
    
    lukM=57,
    
    lukMF=58,
    
    lukPVF=59,
    
    lukPVS=60,
    
    sak=61,
    
    sasX=62,
    
    scn=63,
    
    sea=64,
    
    seb=65,
    
    sec=66,
    
    sed=67,
    
    see=68,
    
    seg=69,
    
    seh=70,
    
    sei=71,
    
    sej=72,
    
    selR=73,
    
    sep=74,
    
    seu=75,
    
    tsst1=76,
    
    unspecified_gpg = 77
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 77    //ignore unspecified_gpg
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