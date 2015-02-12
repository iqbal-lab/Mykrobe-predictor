/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/

/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  base_species.h
*/
#include "dB_graph.h"

#ifdef STAPH

typedef enum 
 {
   Staphylococcus = 0  
} Complex ;

#define NUM_COMPLEX 2

 typedef enum 
{
Saureus=0,
Sepidermidis = 1,
Shaemolyticus = 2,
Sother=3
} Species ;

#define NUM_SPECIES 4


#define NUM_LINEAGE 0

 #endif

#ifdef TB


typedef enum 
 {
   MTBC = 0,
   NTM = 1
  } Complex ;

#define NUM_COMPLEX 2

typedef enum 
 {
  abscessus = 0,
  africanum = 1,
  aromaticivorans = 2,
  avium = 3,
  bovis = 4,
  branderi = 5,
  caprae = 6,
  chelonae = 7,
  chlorophenolicum = 8,
  chubuense = 9,
  colombiense = 10,
  crocinum = 11,
  flavescens = 12,
  fluoranthenivorans = 13,
  fortuitum = 14,
  gilvum = 15,
  gordonae = 16,
  hodleri = 17,
  interjectum = 18,
  intracellulare = 19,
  kansasii = 20,
  lentiflavum = 21,
  leprae = 22,
  malmoense = 23,
  marinum = 24,
  mucogenicum = 25, 
  pallens = 26, 
  peregrinum = 27, 
  phage = 28, 
  pyrenivorans = 29, 
  rufum = 30,
  rutilum = 31,
  scrofulaceum = 32,
  senegalense = 33,
  smegmatis = 34,
  sphagni = 35,
  szulgai = 36,
  triplex = 37,
  tuberculosis = 38,
  tusciae = 39,
  ulcerans = 40,
  vaccae = 41,
  xenopi = 42
  } Species ;

#define NUM_SPECIES 43

  typedef enum 
 {
  beijing = 0,
  lineage1 = 1,
  lineage2 = 2, 
  lineage3 = 3, 
  lineage4 = 4, 
  lineage5 = 5
  } Lineage ;
 #define NUM_LINEAGES 6

 #endif
typedef struct
{
  boolean present[NUM_SPECIES];
  int percentage_coverage[NUM_SPECIES];
  int median_coverage[NUM_SPECIES];
  int percentage_threshold[NUM_SPECIES];  
  int num_panels_present;
  int NUM_PANELS;
} CovgInfo;



typedef struct
{
  CovgInfo* complex_covg_info;
  CovgInfo* species_covg_info;
  CovgInfo* lineage_covg_info;
} SpeciesInfo;

void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                            StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first, int ignore_last , int NUM_PANELS);

int get_ith_coverage_panel(CovgInfo* covg_info, int i);
int get_ith_present_panel(CovgInfo* covg_info, int i);


boolean is_percentage_coverage_above_threshold(int per_cov,int threshold);


CovgInfo* alloc_and_init_covg_info();
CovgInfo* get_coverage_info(dBGraph *db_graph,
                          StrBuf** file_paths,
                          int max_branch_len,
                          int NUM_PANELS,
                          int ignore_first,
                          int ignore_last,
                          void (*load_thresholds)()
                          );

int max(int int1,int int2);
int get_best_hit(CovgInfo* covg_info,boolean* mask);

boolean panels_are_present(CovgInfo* covg_info ,  boolean* mask);

void print_json_phylogenetics(SpeciesInfo* species_info);
void print_json_complex(SpeciesInfo* species_info);
void print_json_species(SpeciesInfo* species_info);
void print_json_lineage(SpeciesInfo* species_info);
