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
    
    mphC=32,
    
    sat4=33,
    
    cat=34,
    
    vanA=35,
    
    vanC=36,
    
    vanB=37,
    
    fusB=38,
    
    fusC=39,
    
    str=40,
    
    qacCsmr=41,
    
    arcA=42,
    
    arcB=43,
    
    arcC=44,
    
    arcD=45,
    
    ccrA=46,
    
    ccrB=47,
    
    ccrCa=48,
    
    ccrCb=49,
    
    ccrCc=50,
    
    chp=51,
    
    eta=52,
    
    etb=53,
    
    etd=54,
    
    luk=55,
    
    lukM=56,
    
    lukMF=57,
    
    lukPVF=58,
    
    lukPVS=59,
    
    sak=60,
    
    sasX=61,
    
    scn=62,
    
    sea=63,
    
    seb=64,
    
    sec=65,
    
    sed=66,
    
    see=67,
    
    seg=68,
    
    seh=69,
    
    sei=70,
    
    sej=71,
    
    selR=72,
    
    sep=73,
    
    seu=74,
    
    tsst1=75,
    
    unspecified_gpg = 76
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 76    //ignore unspecified_gpg
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