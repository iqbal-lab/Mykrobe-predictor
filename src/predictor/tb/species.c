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
#include "species.h"
#include "gene_presence.h"
#include "genotyping_known.h"


void map_complex_enum_to_str(Myc_complex sp, StrBuf* sbuf)
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


void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf)
{
  if (sp==abscessus)
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

void map_lineage_enum_to_str(Myc_lineage sp, StrBuf* sbuf)
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



int get_ith_coverage_panel(CovgInfo* covg_info, int i)
{
  int covg=0;
  int j=0;
  int panel_index=0;
  if (i > covg_info->num_panels_present -1 ){
    die("We only have %i species, we can't find %ith coverage \n",  covg_info->num_panels_present,i+1);
  }
  else
  {
    for (j=0; j<covg_info->NUM_PANELS; j++)
    {
      if (covg_info->present[j])
      {
        if (panel_index == i)
        {
            // We at the required species
            covg = covg_info->median_coverage[j];
        }
        panel_index = panel_index + 1;
      }
    }    
  }
  return covg;

}
int get_ith_present_panel(CovgInfo* covg_info, int i){
  int  panel=10000000;
  int j=0; 
  int panel_index = 0;  
  if (i > covg_info->num_panels_present -1 ){
    die("We only have %i panel, we can't find %ith present \n",  covg_info->num_panels_present,i+1);
  }
  else
  {
    for (j=0; j<covg_info->NUM_PANELS; j++)
    {
      if (covg_info->present[j])
      {
        if (panel_index == i)
        {
            // We at the required panel
            panel = j;
        }
        panel_index = panel_index + 1;
      }
    }    
  }
  return (panel);
}

char* get_char_name_of_species_enum(Myc_species species){
  StrBuf* species_name = strbuf_new(); 
  map_species_enum_to_str(species, species_name);
  return species_name->buff;
}

char* get_char_name_of_lineage_enum(Myc_lineage lineage){
  StrBuf* lineage_name = strbuf_new(); 
  map_lineage_enum_to_str(lineage, lineage_name);
  return lineage_name->buff;
}

char* get_ith_complex_name(CovgInfo* covg_info, int i)
{
  Myc_complex complex;
  StrBuf* complex_name = strbuf_new(); 
  complex = get_ith_present_panel( covg_info, i);
  map_complex_enum_to_str(complex, complex_name);
  return complex_name->buff;
}
char* get_ith_species_name(CovgInfo* covg_info, int i)
{
  Myc_species species;
  species = get_ith_present_panel( covg_info, i);
  return get_char_name_of_species_enum(species);
}

char* get_ith_lineage_name(CovgInfo* covg_info, int i)
{
  Myc_lineage lineage;
  lineage = get_ith_present_panel( covg_info, i);
  return get_char_name_of_lineage_enum(lineage);
}

void get_coverage_on_panels(int* percentage_coverage,int* median_coverage,
                            StrBuf** panel_file_paths,
                            int max_branch_len, dBGraph *db_graph,
                            int ignore_first, int ignore_last , int NUM_PANELS)

{
  int i;
  FILE* fp;
  AlleleInfo* ai = alloc_allele_info();
  int number_of_reads = 0;
  int number_of_covered_reads =0; 
  int num_kmers=0;
  Covg tot_pos_kmers;
  Covg tot_kmers;
  Covg med;  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL)
  {
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_branch_len,LINE_MAX);
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
  {
    die("Failed to malloc kmer sliding window");
  }
  CovgArray* working_ca = alloc_and_init_covg_array(max_branch_len);
  dBNode** array_nodes = (dBNode**) malloc(sizeof(dBNode*)*max_branch_len);
  Orientation* array_or =(Orientation*)  malloc(sizeof(Orientation)*max_branch_len);
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_branch_len-db_graph->kmer_size+1));
  if (kmer_window->kmer==NULL)
  {
    die("Failed to malloc kmer_window->kmer");
  }
  kmer_window->nkmers=0;
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false)
    {
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    return ret;
  }
  if ( (array_nodes==NULL) || (array_or==NULL))
  {
    die("Cannot alloc array of nodes or of orientations");
  }
  for (i = 0; i < NUM_PANELS; i++)
  {
    // Reset counters
    tot_pos_kmers = 0; 
    tot_kmers = 0; 
    med = 0;
    number_of_covered_reads = 0;
      fp = fopen(panel_file_paths[i]->buff, "r");
      if (fp==NULL)
      {
        die("Cannot open this file - %s", panel_file_paths[i]->buff);
      } 
      // while the entry is valid iterate through the fasta file
      do 
      {
        num_kmers= get_next_single_allele_info(fp, db_graph, ai,
                                               true,
                                               seq, kmer_window,
                                               &file_reader_fasta,
                                               array_nodes, array_or, 
                                               working_ca, max_branch_len,
                                               ignore_first, ignore_last);

        number_of_reads = number_of_reads + 1;
        int pos_kmers = num_kmers * (double) (ai->percent_nonzero)/100;
        if (ai->percent_nonzero > 0 ){
          number_of_covered_reads = number_of_covered_reads +1; 
        }
        //calculate a running pseudo median, before you update the tots
        if  (tot_kmers+num_kmers>0)
          {
              med += ai->median_covg ; 
              tot_kmers += num_kmers;
              tot_pos_kmers += pos_kmers;
          }
      }while ( num_kmers > 0);

    if ( (number_of_reads > 0) && (tot_kmers>0) )
    {
      percentage_coverage[i] = (int) (100 * tot_pos_kmers) / tot_kmers;
      median_coverage[i] = (double) (med) / number_of_covered_reads ;
    }
        else
    {
      percentage_coverage[i]=0;
      median_coverage[i]=0;
    } 
  fclose(fp);
  }
}

boolean is_percentage_coverage_above_threshold(int per_cov,int threshold)
{
  if (per_cov >= threshold)
  {
    return (true);
  }
  else
  {
    return (false);
  }
}

void find_which_panels_are_present(CovgInfo* covg_info)
{
  int i;
  for (i=0; i<covg_info->NUM_PANELS; i++)
  {
    // printf("panel %d has coverage %d \n",i,covg_info->percentage_coverage[i]);
    covg_info->present[i] = is_percentage_coverage_above_threshold(covg_info->percentage_coverage[i],covg_info->percentage_threshold[i]);
    if (covg_info->present[i])
    {
      covg_info->num_panels_present = covg_info->num_panels_present +1;
    }
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
CovgInfo* alloc_and_init_covg_info(){
  CovgInfo* covg_info=(CovgInfo *)malloc(sizeof(CovgInfo));   
  covg_info->num_panels_present = 0 ;
  int j;
  for(j = 0; j < NUM_SPECIES; j++) {
    covg_info->percentage_coverage[j] = 0;
    covg_info->median_coverage[j] =0;
    covg_info->present[j]=false;
  }    
  return (covg_info);
}
CovgInfo* get_coverage_info(dBGraph *db_graph,
                          StrBuf** file_paths,
                          int max_branch_len,
                          int NUM_PANELS,
                          int ignore_first,
                          int ignore_last,
                          void (*load_thresholds)()
                          ){

  CovgInfo* covg_info = alloc_and_init_covg_info();
  covg_info->NUM_PANELS = NUM_PANELS;  
  load_thresholds(covg_info->percentage_threshold);
  get_coverage_on_panels(covg_info->percentage_coverage,
                          covg_info->median_coverage,
                          file_paths,
                          max_branch_len,db_graph,
                          ignore_first,ignore_last,
                          NUM_PANELS);
  find_which_panels_are_present(covg_info);  

  return (covg_info);
}
void load_all_complex_thresholds(int* thresholds){
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

int max(int int1,int int2){
  if (int1 > int2){
    return (int1);
  }
  else{
    return (int2);
  }

}

void update_complex_presence_and_coverage_from_species(SpeciesInfo* species_info){
  if (MTBC_panels_are_present(species_info)){
    if (! species_info->complex_covg_info->present[MTBC]){
      species_info->complex_covg_info->present[MTBC] = true;
      species_info->complex_covg_info->num_panels_present = species_info->complex_covg_info->num_panels_present + 1;
    }
    Myc_species best_MTBC_species = get_best_MTBC_species(species_info);
    species_info->complex_covg_info->percentage_coverage[MTBC] = max(species_info->complex_covg_info->percentage_coverage[MTBC] , species_info->species_covg_info->percentage_coverage[best_MTBC_species] );
    species_info->complex_covg_info->median_coverage[MTBC] = max(species_info->complex_covg_info->median_coverage[MTBC] , species_info->species_covg_info->median_coverage[best_MTBC_species] );
  }
  if (NTM_panels_are_present(species_info)){
    if (! species_info->complex_covg_info->present[NTM]){
      species_info->complex_covg_info->present[NTM] = true;
      species_info->complex_covg_info->num_panels_present = species_info->complex_covg_info->num_panels_present + 1;
    }
    Myc_species best_NTM_species = get_best_NTM_species(species_info);
    species_info->complex_covg_info->percentage_coverage[NTM] = max(species_info->complex_covg_info->percentage_coverage[NTM] , species_info->species_covg_info->percentage_coverage[best_NTM_species] );
    species_info->complex_covg_info->median_coverage[NTM] = max(species_info->complex_covg_info->median_coverage[NTM] , species_info->species_covg_info->median_coverage[best_NTM_species] );
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


  CovgInfo* complex_covg_info = get_coverage_info(db_graph,
                                                  mtbc_and_ntm_file_paths,
                                                  max_branch_len,NUM_COMPLEX,
                                                  ignore_first,ignore_last,
                                                  load_all_complex_thresholds);
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
  species_info->complex_covg_info = complex_covg_info;
  species_info->species_covg_info = species_covg_info;
  species_info->lineage_covg_info = lineage_covg_info;

  update_complex_presence_and_coverage_from_species(species_info);

  return species_info;
}

// boolean sample_is_mixed(boolean NTM_is_present,boolean MTBC_is_present)

// {
//   boolean is_mixed = NTM_is_present && MTBC_is_present;
//   return (is_mixed);
// }

// boolean sample_is_MTBC(boolean NTM_is_present,boolean MTBC_is_present)
// {
//   boolean is_MTBC = MTBC_is_present && !NTM_is_present;
//   return (is_MTBC);
// }

// boolean sample_is_NTM(boolean NTM_is_present,boolean MTBC_is_present)
// {
//   boolean is_NTM = !MTBC_is_present && NTM_is_present;
//   return (is_NTM);
// }

boolean is_NTM_present(SpeciesInfo* species_info)
{
  // Is combined NTM panel present OR any of the NTM species panels
  if (species_info->complex_covg_info->present[NTM])
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

  if (species_info->complex_covg_info->present[MTBC]){
    return (true);
  }
  else{
    return (false);
  }
}


// SampleType get_sample_type(SpeciesInfo* species_info)
// {
//   boolean NTM_is_present = is_NTM_present(species_info);
//   boolean MTBC_is_present = is_MTBC_present(species_info);
//   SampleType sample_type;
//   if ( sample_is_mixed(NTM_is_present,MTBC_is_present) )
//     {
//       sample_type=MixedTB;
//     }
//   else{
//     if ( sample_is_MTBC(NTM_is_present,MTBC_is_present) ){
//       sample_type = PureMTBC;
//     }
//     else if ( sample_is_NTM(NTM_is_present,MTBC_is_present) )
//     {
//       sample_type = PureNTM;
//     }
//     else
//     {
//       sample_type = NonTB;
//     }    

//   }
//   return sample_type;  
// }

void print_json_indiv_phylo(CovgInfo* covg_info,
                           char* (*get_ith_name)(CovgInfo*, int)){
    int i;
    boolean last = false;
    for (i=0; i < covg_info->num_panels_present; i++)
    {
      if (i == covg_info->num_panels_present-1){
        last = true;
      }
      print_json_called_variant_item( (*get_ith_name)(covg_info,i), get_ith_coverage_panel(covg_info,i), last);
    }     
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
      print_json_called_variant_item( "Non Mycobacterium", -1, true);
    }
    print_json_complex_end();  
}

Myc_species get_best_hit(CovgInfo* covg_info,boolean* mask)
{
  int i;
  int best_perc_cov_so_far=0;
  int best_median_cov_so_far=0;
  int curr=-1;

  for (i=0; i<covg_info->NUM_PANELS; i++)
  {
    if (mask[i]){
      if (covg_info->percentage_coverage[i] >= best_perc_cov_so_far)
      {
         best_perc_cov_so_far = covg_info->percentage_coverage[i];
        // Only update if the median coverage has also improved
        if (covg_info->median_coverage[i] > best_median_cov_so_far){
          best_median_cov_so_far = covg_info->median_coverage[i];
          curr=i;
        }          
      }      
    }
  }
  return (Myc_lineage) curr;
}

boolean* create_mask(boolean default_value)
{
  boolean* mask = malloc(NUM_SPECIES * sizeof(boolean));
  int j;
  for(j = 0; j < NUM_SPECIES; j++) {
    mask[j] = default_value;
  }   
  return (mask);
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

Myc_species get_best_MTBC_species(SpeciesInfo* species_info ){
  boolean* mask = create_MTBC_mask();
  Myc_species species = get_best_hit(species_info->species_covg_info,mask);
  return (species);
}

Myc_species get_best_NTM_species(SpeciesInfo* species_info ){
  boolean* mask = create_NTM_mask();
  Myc_species species = get_best_hit(species_info->species_covg_info,mask);
  return (species);
}

Myc_lineage get_best_lineage(SpeciesInfo* species_info ){
  boolean* mask= create_mask(true);
  Myc_lineage lineage = get_best_hit(species_info->lineage_covg_info,mask);
  return (lineage);
}

boolean panels_are_present(CovgInfo* covg_info ,  boolean* mask){
  boolean panels_are_present = false;
  int i;
  for (i=0; i<covg_info->NUM_PANELS; i++)
  {
    if (mask[i]){
      if (covg_info->present[i]){
        panels_are_present = true;
      }
    }
  }
  return (panels_are_present);
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
    print_json_called_variant_item( "Unknown MTBC Species", -1 , false);
  }
  else
  {  
    Myc_species MTBC_species = get_best_MTBC_species(species_info);
    print_json_called_variant_item( get_char_name_of_species_enum(MTBC_species), species_info->species_covg_info->median_coverage[MTBC_species], false);
  }
  if (no_NTM_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown NTM Species", -1 , true);
  }
  else
  {  
    Myc_species NTM_species = get_best_NTM_species(species_info);  
    print_json_called_variant_item( get_char_name_of_species_enum(NTM_species), species_info->species_covg_info->median_coverage[NTM_species], true);
  }
}

void print_json_best_hit_MBTC(SpeciesInfo* species_info){
  if (no_MTBC_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown MTBC Species", -1 , true);
  }
  else{
  Myc_species MTBC_species = get_best_MTBC_species(species_info);
  print_json_called_variant_item( get_char_name_of_species_enum(MTBC_species), species_info->species_covg_info->median_coverage[MTBC_species], true);    
  }
}

void print_json_best_hit_NTM(SpeciesInfo* species_info){
  if (no_NTM_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown NTM Species", -1 , true);
  }
  else
  {  
    Myc_species NTM_species = get_best_NTM_species(species_info);  
    print_json_called_variant_item( get_char_name_of_species_enum(NTM_species), species_info->species_covg_info->median_coverage[NTM_species], true);
  }

}

void print_json_best_hit_lineage(SpeciesInfo* species_info){
  if (no_lineage_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Lineage", -1 , true);
  }
  else
  {    
    Myc_lineage lineage = get_best_lineage(species_info);  
    print_json_called_variant_item( get_char_name_of_lineage_enum(lineage), species_info->lineage_covg_info->median_coverage[lineage], true);
  }
}

void print_json_species(SpeciesInfo* species_info){
    Myc_species MTBC_is_present = is_MTBC_present(species_info);
    Myc_species NTM_is_present = is_NTM_present(species_info);
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
  boolean MTBC_is_present = species_info->complex_covg_info->present[MTBC];
  boolean NTM_is_present = species_info->complex_covg_info->present[NTM];
  return (MTBC_is_present || NTM_is_present);
}

void print_json_phylogenetics(SpeciesInfo* species_info){
    print_json_phylogenetics_start();
    print_json_complex(species_info);
    // if (myco_is_present(species_info)){
    print_json_species(species_info);
      // if (tuberculosis_is_present(species_info)){
    print_json_lineage(species_info);
      // }      
    // }
    print_json_phylogenetics_end();  
}
