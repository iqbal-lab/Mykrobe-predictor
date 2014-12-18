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
    NonTB = 2,
  } SampleType;

typedef enum 
 {
   Mtuberculosis = 0,
   Mafricanum = 1,
   Mbovis = 2 
  } Myc_species ;

#define NUM_SPECIES 3

typedef struct
{
  SampleType sample_type;
  boolean present[NUM_SPECIES];
  int percentage_coverage[NUM_SPECIES];
  int median_coverage[NUM_SPECIES];
  int num_species;

} SpeciesInfo;



typedef enum 
 {
  beijing = 0,
  animal = 1, 
  lineage1 = 2,
  lineage2 = 3, 
  lineage3 = 4, 
  lineage4 = 5, 
  lineage5 = 6, 
  lineage6 = 7, 
  lineage7 = 8
  } Myc_lineage ;
 #define NUM_LINEAGES 9



void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                             StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first,int ignore_last,
                            int NUM_PANELS);
boolean is_percentage_coverage_above_threshold(int per_cov,int threshold);
void are_mtbc_and_ntm_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels);
boolean sample_is_mixed(boolean* sample_type_panels);
boolean sample_is_MTBC(boolean* sample_type_panels);
boolean sample_is_NTM(boolean* sample_type_panels);
SampleType get_sample_type(boolean* sample_type_panels, boolean* lineage_panels);

void load_all_mtbc_and_ntm_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );

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


