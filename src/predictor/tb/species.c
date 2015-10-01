/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
  species.c
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
#include "species.h"
#include "gene_presence.h"
#include "genotyping_known.h"


  void map_complex_enum_to_str(Complex sp, StrBuf* sbuf)
{
  
     if(sp==Mtbc){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "MTBC");
      }
  
     else if(sp==Ntm){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "NTM");
      }
  
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}
  void load_all_complex_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
	
  panel_file_paths[Mtbc] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mtbc], "data/tb/phylo/complex/MTBC.fa" );
  
  panel_file_paths[Ntm] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ntm], "data/tb/phylo/complex/NTM.fa" );
  
}
  void complex_threshold(int* thresholds){
	
	  thresholds[Mtbc] = 50;
	  thresholds[Mtbc] = 50;
	
	  thresholds[Ntm] = 50;
	  thresholds[Ntm] = 50;
	
}
  
void print_json_complex(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->complex_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_complex_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_complex_name);
    }
    else
    {
      print_json_called_variant_item( "Unknown complex", -1, true);
    }
    print_json_phylo_group_end();  
}
  char* get_ith_complex_name(CovgInfo* covg_info, int i)
{
  Complex complex;
  StrBuf* complex_name = strbuf_new(); 
  complex = get_ith_present_panel( covg_info, i);
  map_complex_enum_to_str(complex, complex_name);
  return complex_name->buff;
}

  void map_lineage_enum_to_str(Lineage sp, StrBuf* sbuf)
{
  
     if(sp==Lineage1){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lineage1");
      }
  
     else if(sp==Lineage2){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lineage2");
      }
  
     else if(sp==Lineage3){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lineage3");
      }
  
     else if(sp==Lineage4){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lineage4");
      }
  
     else if(sp==Lineage5){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lineage5");
      }
  
     else if(sp==Lineage6){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lineage6");
      }
  
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}
  void load_all_lineage_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
	
  panel_file_paths[Lineage1] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lineage1], "data/tb/phylo/lineage/lineage_1.fa" );
  
  panel_file_paths[Lineage2] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lineage2], "data/tb/phylo/lineage/lineage_2.fa" );
  
  panel_file_paths[Lineage3] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lineage3], "data/tb/phylo/lineage/lineage_3.fa" );
  
  panel_file_paths[Lineage4] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lineage4], "data/tb/phylo/lineage/lineage_4.fa" );
  
  panel_file_paths[Lineage5] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lineage5], "data/tb/phylo/lineage/lineage_5.fa" );
  
  panel_file_paths[Lineage6] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lineage6], "data/tb/phylo/lineage/lineage_6.fa" );
  
}
  void lineage_threshold(int* thresholds){
	
	  thresholds[Lineage1] = 50;
	  thresholds[Lineage1] = 50;
	
	  thresholds[Lineage2] = 50;
	  thresholds[Lineage2] = 50;
	
	  thresholds[Lineage3] = 50;
	  thresholds[Lineage3] = 50;
	
	  thresholds[Lineage4] = 50;
	  thresholds[Lineage4] = 50;
	
	  thresholds[Lineage5] = 50;
	  thresholds[Lineage5] = 50;
	
	  thresholds[Lineage6] = 50;
	  thresholds[Lineage6] = 50;
	
}
  
void print_json_lineage(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->lineage_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_lineage_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_lineage_name);
    }
    else
    {
      print_json_called_variant_item( "Unknown lineage", -1, true);
    }
    print_json_phylo_group_end();  
}
  char* get_ith_lineage_name(CovgInfo* covg_info, int i)
{
  Lineage lineage;
  StrBuf* lineage_name = strbuf_new(); 
  lineage = get_ith_present_panel( covg_info, i);
  map_lineage_enum_to_str(lineage, lineage_name);
  return lineage_name->buff;
}

  void map_species_enum_to_str(Species sp, StrBuf* sbuf)
{
  
     if(sp==Abscessus){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "abscessus");
      }
  
     else if(sp==Africanum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "africanum");
      }
  
     else if(sp==Animallineage){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "animallineage");
      }
  
     else if(sp==Aromaticivorans){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "aromaticivorans");
      }
  
     else if(sp==Avium){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "avium");
      }
  
     else if(sp==Beijingsublineage){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "beijingsublineage");
      }
  
     else if(sp==Bovis){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "bovis");
      }
  
     else if(sp==Branderi){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "branderi");
      }
  
     else if(sp==Caprae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "caprae");
      }
  
     else if(sp==Chelonae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "chelonae");
      }
  
     else if(sp==Chlorophenolicum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "chlorophenolicum");
      }
  
     else if(sp==Chubuense){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "chubuense");
      }
  
     else if(sp==Colombiense){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "colombiense");
      }
  
     else if(sp==Crocinum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "crocinum");
      }
  
     else if(sp==Flavescens){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "flavescens");
      }
  
     else if(sp==Fluoranthenivorans){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "fluoranthenivorans");
      }
  
     else if(sp==Fortuitum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "fortuitum");
      }
  
     else if(sp==Gilvum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "gilvum");
      }
  
     else if(sp==Gordonae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "gordonae");
      }
  
     else if(sp==Hodleri){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "hodleri");
      }
  
     else if(sp==Interjectum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "interjectum");
      }
  
     else if(sp==Intracellulare){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "intracellulare");
      }
  
     else if(sp==Kansasii){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "kansasii");
      }
  
     else if(sp==Lentiflavum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "lentiflavum");
      }
  
     else if(sp==Leprae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "leprae");
      }
  
     else if(sp==Malmoense){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "malmoense");
      }
  
     else if(sp==Marinum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "marinum");
      }
  
     else if(sp==Mucogenicum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "mucogenicum");
      }
  
     else if(sp==Pallens){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "pallens");
      }
  
     else if(sp==Peregrinum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "peregrinum");
      }
  
     else if(sp==Phage){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "phage");
      }
  
     else if(sp==Pyrenivorans){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "pyrenivorans");
      }
  
     else if(sp==Rhodesiae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "rhodesiae");
      }
  
     else if(sp==Rufum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "rufum");
      }
  
     else if(sp==Rutilum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "rutilum");
      }
  
     else if(sp==Scrofulaceum){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "scrofulaceum");
      }
  
     else if(sp==Senegalense){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "senegalense");
      }
  
     else if(sp==Smegmatis){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "smegmatis");
      }
  
     else if(sp==Sphagni){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "sphagni");
      }
  
     else if(sp==Szulgai){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "szulgai");
      }
  
     else if(sp==Triplex){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "triplex");
      }
  
     else if(sp==Tuberculosis){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "tuberculosis");
      }
  
     else if(sp==Tusciae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "tusciae");
      }
  
     else if(sp==Ulcerans){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "ulcerans");
      }
  
     else if(sp==Vaccae){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "vaccae");
      }
  
     else if(sp==Xenopi){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "xenopi");
      }
  
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}
  void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
	
  panel_file_paths[Abscessus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Abscessus], "data/tb/phylo/species/abscessus.fa" );
  
  panel_file_paths[Africanum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Africanum], "data/tb/phylo/species/africanum.fa" );
  
  panel_file_paths[Animallineage] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Animallineage], "data/tb/phylo/species/animal_lineage.fa" );
  
  panel_file_paths[Aromaticivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aromaticivorans], "data/tb/phylo/species/aromaticivorans.fa" );
  
  panel_file_paths[Avium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Avium], "data/tb/phylo/species/avium.fa" );
  
  panel_file_paths[Beijingsublineage] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Beijingsublineage], "data/tb/phylo/species/beijing_sublineage.fa" );
  
  panel_file_paths[Bovis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bovis], "data/tb/phylo/species/bovis.fa" );
  
  panel_file_paths[Branderi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Branderi], "data/tb/phylo/species/branderi.fa" );
  
  panel_file_paths[Caprae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caprae], "data/tb/phylo/species/caprae.fa" );
  
  panel_file_paths[Chelonae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chelonae], "data/tb/phylo/species/chelonae.fa" );
  
  panel_file_paths[Chlorophenolicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorophenolicum], "data/tb/phylo/species/chlorophenolicum.fa" );
  
  panel_file_paths[Chubuense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chubuense], "data/tb/phylo/species/chubuense.fa" );
  
  panel_file_paths[Colombiense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Colombiense], "data/tb/phylo/species/colombiense.fa" );
  
  panel_file_paths[Crocinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Crocinum], "data/tb/phylo/species/crocinum.fa" );
  
  panel_file_paths[Flavescens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Flavescens], "data/tb/phylo/species/flavescens.fa" );
  
  panel_file_paths[Fluoranthenivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fluoranthenivorans], "data/tb/phylo/species/fluoranthenivorans.fa" );
  
  panel_file_paths[Fortuitum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fortuitum], "data/tb/phylo/species/fortuitum.fa" );
  
  panel_file_paths[Gilvum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Gilvum], "data/tb/phylo/species/gilvum.fa" );
  
  panel_file_paths[Gordonae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Gordonae], "data/tb/phylo/species/gordonae.fa" );
  
  panel_file_paths[Hodleri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Hodleri], "data/tb/phylo/species/hodleri.fa" );
  
  panel_file_paths[Interjectum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Interjectum], "data/tb/phylo/species/interjectum.fa" );
  
  panel_file_paths[Intracellulare] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Intracellulare], "data/tb/phylo/species/intracellulare.fa" );
  
  panel_file_paths[Kansasii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Kansasii], "data/tb/phylo/species/kansasii.fa" );
  
  panel_file_paths[Lentiflavum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lentiflavum], "data/tb/phylo/species/lentiflavum.fa" );
  
  panel_file_paths[Leprae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leprae], "data/tb/phylo/species/leprae.fa" );
  
  panel_file_paths[Malmoense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Malmoense], "data/tb/phylo/species/malmoense.fa" );
  
  panel_file_paths[Marinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Marinum], "data/tb/phylo/species/marinum.fa" );
  
  panel_file_paths[Mucogenicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mucogenicum], "data/tb/phylo/species/mucogenicum.fa" );
  
  panel_file_paths[Pallens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pallens], "data/tb/phylo/species/pallens.fa" );
  
  panel_file_paths[Peregrinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Peregrinum], "data/tb/phylo/species/peregrinum.fa" );
  
  panel_file_paths[Phage] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Phage], "data/tb/phylo/species/phage.fa" );
  
  panel_file_paths[Pyrenivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrenivorans], "data/tb/phylo/species/pyrenivorans.fa" );
  
  panel_file_paths[Rhodesiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhodesiae], "data/tb/phylo/species/rhodesiae.fa" );
  
  panel_file_paths[Rufum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rufum], "data/tb/phylo/species/rufum.fa" );
  
  panel_file_paths[Rutilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rutilum], "data/tb/phylo/species/rutilum.fa" );
  
  panel_file_paths[Scrofulaceum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Scrofulaceum], "data/tb/phylo/species/scrofulaceum.fa" );
  
  panel_file_paths[Senegalense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Senegalense], "data/tb/phylo/species/senegalense.fa" );
  
  panel_file_paths[Smegmatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Smegmatis], "data/tb/phylo/species/smegmatis.fa" );
  
  panel_file_paths[Sphagni] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sphagni], "data/tb/phylo/species/sphagni.fa" );
  
  panel_file_paths[Szulgai] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Szulgai], "data/tb/phylo/species/szulgai.fa" );
  
  panel_file_paths[Triplex] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Triplex], "data/tb/phylo/species/triplex.fa" );
  
  panel_file_paths[Tuberculosis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Tuberculosis], "data/tb/phylo/species/tuberculosis.fa" );
  
  panel_file_paths[Tusciae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Tusciae], "data/tb/phylo/species/tusciae.fa" );
  
  panel_file_paths[Ulcerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ulcerans], "data/tb/phylo/species/ulcerans.fa" );
  
  panel_file_paths[Vaccae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vaccae], "data/tb/phylo/species/vaccae.fa" );
  
  panel_file_paths[Xenopi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xenopi], "data/tb/phylo/species/xenopi.fa" );
  
}
  void species_threshold(int* thresholds){
	
	  thresholds[Abscessus] = 50;
	  thresholds[Abscessus] = 50;
	
	  thresholds[Africanum] = 50;
	  thresholds[Africanum] = 50;
	
	  thresholds[Animallineage] = 50;
	  thresholds[Animallineage] = 50;
	
	  thresholds[Aromaticivorans] = 50;
	  thresholds[Aromaticivorans] = 50;
	
	  thresholds[Avium] = 50;
	  thresholds[Avium] = 50;
	
	  thresholds[Beijingsublineage] = 50;
	  thresholds[Beijingsublineage] = 50;
	
	  thresholds[Bovis] = 50;
	  thresholds[Bovis] = 50;
	
	  thresholds[Branderi] = 50;
	  thresholds[Branderi] = 50;
	
	  thresholds[Caprae] = 50;
	  thresholds[Caprae] = 50;
	
	  thresholds[Chelonae] = 50;
	  thresholds[Chelonae] = 50;
	
	  thresholds[Chlorophenolicum] = 50;
	  thresholds[Chlorophenolicum] = 50;
	
	  thresholds[Chubuense] = 50;
	  thresholds[Chubuense] = 50;
	
	  thresholds[Colombiense] = 50;
	  thresholds[Colombiense] = 50;
	
	  thresholds[Crocinum] = 50;
	  thresholds[Crocinum] = 50;
	
	  thresholds[Flavescens] = 50;
	  thresholds[Flavescens] = 50;
	
	  thresholds[Fluoranthenivorans] = 50;
	  thresholds[Fluoranthenivorans] = 50;
	
	  thresholds[Fortuitum] = 50;
	  thresholds[Fortuitum] = 50;
	
	  thresholds[Gilvum] = 50;
	  thresholds[Gilvum] = 50;
	
	  thresholds[Gordonae] = 50;
	  thresholds[Gordonae] = 50;
	
	  thresholds[Hodleri] = 50;
	  thresholds[Hodleri] = 50;
	
	  thresholds[Interjectum] = 50;
	  thresholds[Interjectum] = 50;
	
	  thresholds[Intracellulare] = 50;
	  thresholds[Intracellulare] = 50;
	
	  thresholds[Kansasii] = 50;
	  thresholds[Kansasii] = 50;
	
	  thresholds[Lentiflavum] = 50;
	  thresholds[Lentiflavum] = 50;
	
	  thresholds[Leprae] = 50;
	  thresholds[Leprae] = 50;
	
	  thresholds[Malmoense] = 50;
	  thresholds[Malmoense] = 50;
	
	  thresholds[Marinum] = 50;
	  thresholds[Marinum] = 50;
	
	  thresholds[Mucogenicum] = 50;
	  thresholds[Mucogenicum] = 50;
	
	  thresholds[Pallens] = 50;
	  thresholds[Pallens] = 50;
	
	  thresholds[Peregrinum] = 50;
	  thresholds[Peregrinum] = 50;
	
	  thresholds[Phage] = 50;
	  thresholds[Phage] = 50;
	
	  thresholds[Pyrenivorans] = 50;
	  thresholds[Pyrenivorans] = 50;
	
	  thresholds[Rhodesiae] = 50;
	  thresholds[Rhodesiae] = 50;
	
	  thresholds[Rufum] = 50;
	  thresholds[Rufum] = 50;
	
	  thresholds[Rutilum] = 50;
	  thresholds[Rutilum] = 50;
	
	  thresholds[Scrofulaceum] = 50;
	  thresholds[Scrofulaceum] = 50;
	
	  thresholds[Senegalense] = 50;
	  thresholds[Senegalense] = 50;
	
	  thresholds[Smegmatis] = 50;
	  thresholds[Smegmatis] = 50;
	
	  thresholds[Sphagni] = 50;
	  thresholds[Sphagni] = 50;
	
	  thresholds[Szulgai] = 50;
	  thresholds[Szulgai] = 50;
	
	  thresholds[Triplex] = 50;
	  thresholds[Triplex] = 50;
	
	  thresholds[Tuberculosis] = 50;
	  thresholds[Tuberculosis] = 50;
	
	  thresholds[Tusciae] = 50;
	  thresholds[Tusciae] = 50;
	
	  thresholds[Ulcerans] = 50;
	  thresholds[Ulcerans] = 50;
	
	  thresholds[Vaccae] = 50;
	  thresholds[Vaccae] = 50;
	
	  thresholds[Xenopi] = 50;
	  thresholds[Xenopi] = 50;
	
}
  
void print_json_species(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->species_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_species_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_species_name);
    }
    else
    {
      print_json_called_variant_item( "Unknown species", -1, true);
    }
    print_json_phylo_group_end();  
}
  char* get_ith_species_name(CovgInfo* covg_info, int i)
{
  Species species;
  StrBuf* species_name = strbuf_new(); 
  species = get_ith_present_panel( covg_info, i);
  map_species_enum_to_str(species, species_name);
  return species_name->buff;
}




SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{

  char* get_ith_complex_name(CovgInfo* covg_info, int i);
  StrBuf* complex_file_paths[NUM_Complex];
  load_all_complex_file_paths(complex_file_paths,install_dir);
  CovgInfo* complex_covg_info = get_coverage_info(db_graph,
                                                  complex_file_paths,
                                                  max_branch_len,NUM_Complex,
                                                  ignore_first,ignore_last,
                                                  complex_threshold);  


  char* get_ith_lineage_name(CovgInfo* covg_info, int i);
  StrBuf* lineage_file_paths[NUM_Lineage];
  load_all_lineage_file_paths(lineage_file_paths,install_dir);
  CovgInfo* lineage_covg_info = get_coverage_info(db_graph,
                                                  lineage_file_paths,
                                                  max_branch_len,NUM_Lineage,
                                                  ignore_first,ignore_last,
                                                  lineage_threshold);  


  char* get_ith_species_name(CovgInfo* covg_info, int i);
  StrBuf* species_file_paths[NUM_Species];
  load_all_species_file_paths(species_file_paths,install_dir);
  CovgInfo* species_covg_info = get_coverage_info(db_graph,
                                                  species_file_paths,
                                                  max_branch_len,NUM_Species,
                                                  ignore_first,ignore_last,
                                                  species_threshold);  




  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
    
  species_info->complex_covg_info = complex_covg_info;
    
  species_info->lineage_covg_info = lineage_covg_info;
    
  species_info->species_covg_info = species_covg_info;
  

  // update_phylo_group_presence_and_coverage_from_species(species_info);

  return species_info;
}




void print_json_complex_start()
{
  printf("\t\t\"complex\": {\n");
}


void print_json_lineage_start()
{
  printf("\t\t\"lineage\": {\n");
}


void print_json_species_start()
{
  printf("\t\t\"species\": {\n");
}


void print_json_phylogenetics(SpeciesInfo* species_info){
    print_json_phylogenetics_start();

  
    print_json_complex(species_info);
  
    print_json_lineage(species_info);
  
    print_json_species(species_info);
  
    print_json_phylogenetics_end();  
}

