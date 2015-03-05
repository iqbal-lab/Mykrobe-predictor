/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.c
*/
// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <string_buffer.h>

#include "build.h" 
#include "maths.h" 
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "dB_graph.h"
#include "gene_presence.h"
#include "genotyping_known.h"
 #include "species.h"


void map_phylo_group_enum_to_str(PhyloGroup sp, StrBuf* sbuf)
{
  if(sp==MTBC){
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "MTBC");
  }
   else if(sp== NTM){
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "NTM");
   }
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value\n");
    }
}


void map_species_enum_to_str(Species sp, StrBuf* sbuf)
{
  if (sp==unknown)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Unknown Species");
    }  
  else if (sp==abscessus)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. abscessus");
    }
  else if (sp==africanum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. africanum");
  }
  else if (sp==aromaticivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. aromaticivorans");
  }
  else if (sp==avium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. avium");
  }
  else if (sp==bovis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. bovis");
  }
  else if (sp==branderi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. branderi");
  }
  else if (sp==caprae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. caprae");
  }
  else if (sp==chelonae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. chelonae");
  }
  else if (sp==chlorophenolicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. chlorophenolicum");
  }
  else if (sp==chubuense)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. chubuense");
  }
  else if (sp==colombiense)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. colombiense");
  }
  else if (sp==crocinum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. crocinum");
  }
  else if (sp==flavescens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. flavescens");
  }
  else if (sp==fluoranthenivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. fluoranthenivorans");
  }
  else if (sp==fortuitum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. fortuitum");
  }
  else if (sp==gilvum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. gilvum");
  }
  else if (sp==gordonae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. gordonae");
  }
  else if (sp==hodleri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. hodleri");
  }
  else if (sp==interjectum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. interjectum");
  }
  else if (sp==intracellulare)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. intracellulare");
  }
  else if (sp==kansasii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. kansasii");
  }
  else if (sp==lentiflavum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. lentiflavum");
  }
  else if (sp==leprae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. leprae");
  }
  else if (sp==malmoense)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. malmoense");
  }
  else if (sp==marinum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. marinum");
  }
  else if (sp==mucogenicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. mucogenicum");
  }
  else if (sp==pallens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. pallens");
  }
  else if (sp==peregrinum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. peregrinum");
  }
  else if (sp==phage)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. phage");
  }
  else if (sp==pyrenivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. pyrenivorans");
  }
  else if (sp==rufum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. rufum");
  }
  else if (sp==rutilum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. rutilum");
  }
  else if (sp==scrofulaceum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. scrofulaceum");
  }
  else if (sp==senegalense)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. senegalense");
  }
  else if (sp==smegmatis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. smegmatis");
  }
  else if (sp==sphagni)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. sphagni");
  }
  else if (sp==szulgai)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. szulgai");
  }
  else if (sp==triplex)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. triplex");
  }
  else if (sp==tuberculosis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. tuberculosis");
  }
  else if (sp==tusciae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. tusciae");
  }
  else if (sp==ulcerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. ulcerans");
  }
  else if (sp==vaccae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. vaccae");
  }
  else if (sp==xenopi    )
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf,"M. xenopi");
  }
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value\n");
    }
}

void map_lineage_enum_to_str(Lineage sp, StrBuf* sbuf)
{
  if (sp==beijing)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Beijing/East Asia");
    }
  else if (sp==lineage1)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"East Africa / Indian ocean");
    }
  else if (sp==lineage2)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Beijing/East Asia");
    }
  else if (sp==lineage3)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Delhi/Central Asia");
    }
  else if (sp==lineage4)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"European/American");
    }
  else if (sp==lineage5)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"Ethiopian");
    }
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value\n");
    }
}



void load_all_mtbc_and_ntm_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[0] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[0], "data/tb/species/MTBC.fa" );
  panel_file_paths[1] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[1], "data/tb/species/NTM.fa" );
}

void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[beijing] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[beijing], "data/tb/species/beijing_sublineage.fa");
  panel_file_paths[lineage1] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[lineage1], "data/tb/species/lineage_1.fa");
  panel_file_paths[lineage2] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[lineage2], "data/tb/species/lineage_2.fa");
  panel_file_paths[lineage3] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[lineage3], "data/tb/species/lineage_3.fa");
  panel_file_paths[lineage4] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[lineage4], "data/tb/species/lineage_4.fa");
  panel_file_paths[lineage5] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[lineage5], "data/tb/species/lineage_5.fa");
}

void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[abscessus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[abscessus], "data/tb/species/abscessus.fa");
  panel_file_paths[africanum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[africanum], "data/tb/species/africanum.fa");
  panel_file_paths[aromaticivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[aromaticivorans], "data/tb/species/aromaticivorans.fa");
  panel_file_paths[avium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[avium], "data/tb/species/avium.fa");
  panel_file_paths[bovis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[bovis], "data/tb/species/bovis.fa");
  panel_file_paths[branderi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[branderi], "data/tb/species/branderi.fa");
  panel_file_paths[caprae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[caprae], "data/tb/species/caprae.fa");
  panel_file_paths[chelonae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[chelonae], "data/tb/species/chelonae.fa");
  panel_file_paths[chlorophenolicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[chlorophenolicum], "data/tb/species/chlorophenolicum.fa");
  panel_file_paths[chubuense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[chubuense], "data/tb/species/chubuense.fa");
  panel_file_paths[colombiense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[colombiense], "data/tb/species/colombiense.fa");
  panel_file_paths[crocinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[crocinum], "data/tb/species/crocinum.fa");
  panel_file_paths[flavescens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[flavescens], "data/tb/species/flavescens.fa");
  panel_file_paths[fluoranthenivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[fluoranthenivorans], "data/tb/species/fluoranthenivorans.fa");
  panel_file_paths[fortuitum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[fortuitum], "data/tb/species/fortuitum.fa");
  panel_file_paths[gilvum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[gilvum], "data/tb/species/gilvum.fa");
  panel_file_paths[gordonae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[gordonae], "data/tb/species/gordonae.fa");
  panel_file_paths[hodleri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[hodleri], "data/tb/species/hodleri.fa");
  panel_file_paths[interjectum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[interjectum], "data/tb/species/interjectum.fa");
  panel_file_paths[intracellulare] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[intracellulare], "data/tb/species/intracellulare.fa");
  panel_file_paths[kansasii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[kansasii], "data/tb/species/kansasii.fa");
  panel_file_paths[lentiflavum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[lentiflavum], "data/tb/species/lentiflavum.fa");
  panel_file_paths[leprae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[leprae], "data/tb/species/leprae.fa");
  panel_file_paths[malmoense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[malmoense], "data/tb/species/malmoense.fa");
  panel_file_paths[marinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[marinum], "data/tb/species/marinum.fa");
  panel_file_paths[mucogenicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[mucogenicum], "data/tb/species/mucogenicum.fa");
  panel_file_paths[pallens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[pallens], "data/tb/species/pallens.fa");
  panel_file_paths[peregrinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[peregrinum], "data/tb/species/peregrinum.fa");
  panel_file_paths[phage] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[phage], "data/tb/species/phage.fa");
  panel_file_paths[pyrenivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[pyrenivorans], "data/tb/species/pyrenivorans.fa");
  panel_file_paths[rufum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[rufum], "data/tb/species/rufum.fa");
  panel_file_paths[rutilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[rutilum], "data/tb/species/rutilum.fa");
  panel_file_paths[scrofulaceum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[scrofulaceum], "data/tb/species/scrofulaceum.fa");
  panel_file_paths[senegalense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[senegalense], "data/tb/species/senegalense.fa");
  panel_file_paths[smegmatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[smegmatis], "data/tb/species/smegmatis.fa");
  panel_file_paths[sphagni] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[sphagni], "data/tb/species/sphagni.fa");
  panel_file_paths[szulgai] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[szulgai], "data/tb/species/szulgai.fa");
  panel_file_paths[triplex] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[triplex], "data/tb/species/triplex.fa");
  panel_file_paths[tuberculosis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[tuberculosis], "data/tb/species/tuberculosis.fa");
  panel_file_paths[tusciae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[tusciae], "data/tb/species/tusciae.fa");
  panel_file_paths[ulcerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[ulcerans], "data/tb/species/ulcerans.fa");
  panel_file_paths[vaccae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[vaccae], "data/tb/species/vaccae.fa");
  panel_file_paths[xenopi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[xenopi], "data/tb/species/xenopi.fa");
}

void load_all_phylo_group_thresholds(int* thresholds){
  thresholds[MTBC] = 70;
  thresholds[NTM] = 25;
}
void load_all_species_thresholds(int* thresholds){
  int j;
  for(j = 0; j < NUM_SPECIES; j++) {
    thresholds[j] = 30;
  } 
}
void load_all_lineage_thresholds(int* thresholds){
  int j;
  for(j = 0; j < NUM_SPECIES; j++) {
    thresholds[j] = 90;
  } 
}



void update_phylo_group_presence_and_coverage_from_species(SpeciesInfo* species_info){
  if (MTBC_panels_are_present(species_info)){
    if (! species_info->phylo_group_covg_info->present[MTBC]){
      species_info->phylo_group_covg_info->present[MTBC] = true;
      species_info->phylo_group_covg_info->num_panels_present = species_info->phylo_group_covg_info->num_panels_present + 1;
    }
    Species best_MTBC_species = get_best_MTBC_species(species_info);
    species_info->phylo_group_covg_info->percentage_coverage[MTBC] = max(species_info->phylo_group_covg_info->percentage_coverage[MTBC] , species_info->species_covg_info->percentage_coverage[best_MTBC_species] );
    species_info->phylo_group_covg_info->median_coverage[MTBC] = max(species_info->phylo_group_covg_info->median_coverage[MTBC] , species_info->species_covg_info->median_coverage[best_MTBC_species] );
  }
  if (NTM_panels_are_present(species_info)){
    if (! species_info->phylo_group_covg_info->present[NTM]){
      species_info->phylo_group_covg_info->present[NTM] = true;
      species_info->phylo_group_covg_info->num_panels_present = species_info->phylo_group_covg_info->num_panels_present + 1;
    }
    Species best_NTM_species = get_best_NTM_species(species_info);
    species_info->phylo_group_covg_info->percentage_coverage[NTM] = max(species_info->phylo_group_covg_info->percentage_coverage[NTM] , species_info->species_covg_info->percentage_coverage[best_NTM_species] );
    species_info->phylo_group_covg_info->median_coverage[NTM] = max(species_info->phylo_group_covg_info->median_coverage[NTM] , species_info->species_covg_info->median_coverage[best_NTM_species] );
  }
}
SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{
  StrBuf* mtbc_and_ntm_file_paths[NUM_COMPLEX];
  load_all_mtbc_and_ntm_file_paths(mtbc_and_ntm_file_paths,install_dir);
  StrBuf* species_file_paths[NUM_SPECIES];
  load_all_species_file_paths(species_file_paths,install_dir);
  StrBuf* lineage_file_paths[NUM_LINEAGES];
  load_all_lineage_file_paths(lineage_file_paths,install_dir);


  CovgInfo* phylo_group_covg_info = get_coverage_info(db_graph,
                                                  mtbc_and_ntm_file_paths,
                                                  max_branch_len,NUM_COMPLEX,
                                                  ignore_first,ignore_last,
                                                  load_all_phylo_group_thresholds);
  CovgInfo* species_covg_info = get_coverage_info(db_graph,
                                                  species_file_paths,
                                                  max_branch_len,NUM_SPECIES,
                                                  ignore_first,ignore_last,
                                                  load_all_species_thresholds);
  CovgInfo* lineage_covg_info = get_coverage_info(db_graph,
                                                  lineage_file_paths,
                                                  max_branch_len,NUM_LINEAGES,
                                                  ignore_first,ignore_last,
                                                  load_all_lineage_thresholds);



  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
  species_info->phylo_group_covg_info = phylo_group_covg_info;
  species_info->species_covg_info = species_covg_info;
  species_info->lineage_covg_info = lineage_covg_info;

  update_phylo_group_presence_and_coverage_from_species(species_info);

  return species_info;
}


boolean is_NTM_present(SpeciesInfo* species_info)
{
  // Is combined NTM panel present OR any of the NTM species panels
  if (species_info->phylo_group_covg_info->present[NTM])
  {
    return (true);
  }
  else
  {
    return (false);
  }
}

boolean is_MTBC_present(SpeciesInfo* species_info)
{
  // Is combined MTBC panel present OR any of the MTBC species panels
  if (species_info->phylo_group_covg_info->present[MTBC]){
    return (true);
  }
  else{
    return (false);
  }
}

boolean* create_MTBC_mask()
{
  boolean* mask= create_mask(false);
  mask[tuberculosis] = true;
  mask[bovis] = true;
  mask[caprae] = true;
  mask[africanum] = true;
  return (mask);
}

boolean* create_NTM_mask()
{
  boolean* mask= create_mask(true);
  mask[tuberculosis] = false;
  mask[bovis] = false;
  mask[caprae] = false;
  mask[africanum] = false;
  return (mask);
}

Species get_best_MTBC_species(SpeciesInfo* species_info ){
  boolean* mask = create_MTBC_mask();
  int species_enum = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  return (species);
}

Species get_best_NTM_species(SpeciesInfo* species_info ){
  boolean* mask = create_NTM_mask();
  int species_enum  = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  return (species);
}

Lineage get_best_lineage(SpeciesInfo* species_info ){
  boolean* mask= create_mask(true);
  int lineage_enum  = get_best_hit(species_info->lineage_covg_info,mask);
  Lineage lineage = lineage_enum;
  return (lineage);
}

boolean MTBC_panels_are_present(SpeciesInfo* species_info){
  boolean* mask = create_MTBC_mask();
  boolean MTBC_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  return (MTBC_species_panels_are_present);
}
boolean NTM_panels_are_present(SpeciesInfo* species_info){
  boolean* mask = create_NTM_mask();
  boolean NTM_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  return (NTM_species_panels_are_present);  
}

boolean no_MTBC_panels_are_present(SpeciesInfo* species_info){
  return (!MTBC_panels_are_present(species_info));
}
boolean no_NTM_panels_are_present(SpeciesInfo* species_info){
  return (!NTM_panels_are_present(species_info));  
}

boolean no_lineage_panels_are_present(SpeciesInfo* species_info){
  if (species_info->lineage_covg_info->num_panels_present >0 ){
    return (false);
  }
  else{
    return (true);
  }
}

void print_json_best_hit_NTM_and_MBTC(SpeciesInfo* species_info){
  if (no_MTBC_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , false);
  }
  else
  {  
    Species MTBC_species = get_best_MTBC_species(species_info);
    print_json_called_variant_item( get_char_name_of_species_enum(MTBC_species), species_info->species_covg_info->median_coverage[MTBC_species], false);
  }
  if (no_NTM_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else
  {  
    Species NTM_species = get_best_NTM_species(species_info);  
    print_json_called_variant_item( get_char_name_of_species_enum(NTM_species), species_info->species_covg_info->median_coverage[NTM_species], true);
  }
}

void print_json_best_hit_MBTC(SpeciesInfo* species_info){
  if (no_MTBC_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else{
  Species MTBC_species = get_best_MTBC_species(species_info);
  print_json_called_variant_item( get_char_name_of_species_enum(MTBC_species), species_info->species_covg_info->median_coverage[MTBC_species], true);    
  }
}

void print_json_best_hit_NTM(SpeciesInfo* species_info){
  if (no_NTM_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else
  {  
    Species NTM_species = get_best_NTM_species(species_info);  
    print_json_called_variant_item( get_char_name_of_species_enum(NTM_species), species_info->species_covg_info->median_coverage[NTM_species], true);
  }
}

void print_json_best_hit_lineage(SpeciesInfo* species_info){
  if (no_lineage_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Lineage", -1 , true);
  }
  else
  {    
    Lineage lineage = get_best_lineage(species_info);  
    print_json_called_variant_item( get_char_name_of_lineage_enum(lineage), species_info->lineage_covg_info->median_coverage[lineage], true);
  }
}
void print_json_phylo_group(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->phylo_group_covg_info;
    int num_panels_present = covg_info->num_panels_present;
    print_json_phylo_group_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_phylo_group_name);
    }
    else
    {
      print_json_called_variant_item( "Non Mycobacterium", -1, true);
    }
    print_json_phylo_group_end();  
}

void print_json_species(SpeciesInfo* species_info){
    Species MTBC_is_present = is_MTBC_present(species_info);
    Species NTM_is_present = is_NTM_present(species_info);
    print_json_species_start();
    if (MTBC_is_present && NTM_is_present){
      print_json_best_hit_NTM_and_MBTC(species_info);
    }
    else if (MTBC_is_present){
      print_json_best_hit_MBTC(species_info);

    }
    else if (NTM_is_present){
      print_json_best_hit_NTM(species_info);
    }
    else
    {
      print_json_called_variant_item( "Unknown Species", -1, true);
    }    

    print_json_species_end();  
}
void print_json_lineage(SpeciesInfo* species_info){
    print_json_lineage_start();
    if (tuberculosis_is_present(species_info)){
      print_json_best_hit_lineage(species_info);
    }
    else
    {
      print_json_called_variant_item( "N/A", -1, true);
    }    
    print_json_lineage_end(); 
}

boolean tuberculosis_is_present(SpeciesInfo* species_info){
  return (species_info->species_covg_info->present[tuberculosis]);
}
boolean myco_is_present(SpeciesInfo* species_info){
  boolean MTBC_is_present = species_info->phylo_group_covg_info->present[MTBC];
  boolean NTM_is_present = species_info->phylo_group_covg_info->present[NTM];
  return (MTBC_is_present || NTM_is_present);
}


