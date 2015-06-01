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
    
    ermT=19,
    
    ant9Ib=20,
    
    aphA3aph3III=21,
    
    dfrA=22,
    
    vgaALC=23,
    
    dfrC=24,
    
    dfrD=25,
    
    dfrG=26,
    
    mecC=27,
    
    mecA=28,
    
    dfrK=29,
    
    vgaB=30,
    
    vgaA=31,
    
    tetK=32,
    
    tetM=33,
    
    tetL=34,
    
    tetO=35,
    
    msrA=36,
    
    mphC=37,
    
    sat4=38,
    
    cat=39,
    
    vanA=40,
    
    vanC=41,
    
    vanB=42,
    
    fusB=43,
    
    fusC=44,
    
    str=45,
    
    qacCsmr=46,
    
    arcA=47,
    
    arcB=48,
    
    arcC=49,
    
    arcD=50,
    
    ccrA=51,
    
    ccrB=52,
    
    ccrCa=53,
    
    ccrCb=54,
    
    ccrCc=55,
    
    chp=56,
    
    eta=57,
    
    etb=58,
    
    etd=59,
    
    luk=60,
    
    lukM=61,
    
    lukMF=62,
    
    lukPVF=63,
    
    lukPVS=64,
    
    sak=65,
    
    sasX=66,
    
    scn=67,
    
    sea=68,
    
    seb=69,
    
    sec=70,
    
    sed=71,
    
    see=72,
    
    seg=73,
    
    seh=74,
    
    sei=75,
    
    sej=76,
    
    selR=77,
    
    sep=78,
    
    seu=79,
    
    tsst1=80,
    
    unspecified_gpg = 81
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 81    //ignore unspecified_gpg
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