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
 *  species.h
*/
#include "dB_graph.h"

typedef enum
  {
    PureMTBC =0,
    MixedTB=1,
    PureNTM=2,
    NonTB = 2
  } SampleType;

typedef enum 
 {
   MTBC = 0,
   NTM = 1
  } Myc_complex ;

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
  rhodesiae = 30, 
  rufum = 31, 
  rutilum = 32, 
  scrofulaceum = 33, 
  senegalense = 34, 
  smegmatis = 35, 
  sphagni = 36, 
  szulgai = 37, 
  triplex = 38, 
  tuberculosis = 39, 
  tusciae = 40, 
  ulcerans = 41, 
  vaccae = 42, 
  xenopi = 43
  } Myc_species ;

#define NUM_SPECIES 43

  typedef enum 
 {
  beijing = 0,
  lineage1 = 2,
  lineage2 = 3, 
  lineage3 = 4, 
  lineage4 = 5, 
  lineage5 = 6
  } Myc_lineage ;
 #define NUM_LINEAGES 6


typedef struct
{
  SampleType sample_type;
  boolean present[NUM_SPECIES];
  int percentage_coverage[NUM_SPECIES];
  int median_coverage[NUM_SPECIES];
  int num_species;

} SpeciesInfo;







void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                             StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first,int ignore_last,
                            int NUM_PANELS);
boolean is_percentage_coverage_above_threshold(int per_cov,int threshold);
void are_mtbc_and_ntm_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels);
boolean sample_is_mixed(boolean NTM_is_present,boolean MTBC_is_present);
boolean sample_is_MTBC(boolean NTM_is_present,boolean MTBC_is_present);
boolean sample_is_NTM(boolean NTM_is_present,boolean MTBC_is_present);
boolean is_NTM_present(boolean* complex_presence, boolean* species_presence);
boolean is_MTBC_present(boolean* complex_presence, boolean* species_presence);

SampleType get_sample_type(boolean* complex_presence, boolean* species_presence);

void load_all_mtbc_and_ntm_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );


SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);
void find_which_panels_are_present(int* percentage_coverage,boolean* present,
                                    int* num_panels, boolean has_catalayse);

void map_lineage_enum_to_str(Myc_lineage sp, StrBuf* sbuf);
void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf);

Myc_species map_lineage_enum_to_species_enum(Myc_lineage sp);

char* get_ith_lineage_name(SpeciesInfo* species_info, int i);
char* get_pure_lineage_name(SpeciesInfo* species_info);
int get_ith_lineage_coverage(SpeciesInfo* species_info,int i);
int get_pure_lineage_coverage(SpeciesInfo* species_info);


