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

void map_species_enum_to_str(Myc_species sp, StrBuf* sbuf)
{
  if (sp==tuberculosis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. tuberculosis");
    }
  else if (sp==africanum)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. africanum");
    }
    else if (sp==bovis)
    {
      strbuf_reset(sbuf);
      strbuf_append_str(sbuf,"M. bovis");      
    }

  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value - we get %d\n", sp);
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



char* get_ith_lineage_name(SpeciesInfo* species_info, int i)
{
  Myc_lineage species=0;
  StrBuf* pure_lineage_name = strbuf_new(); 
  int j; 
  int species_index = 0;
  if (i > species_info->num_species -1 ){
    die("We only have %i species, we can't find %i \n",  species_info->num_species,i+1);
  }
  else
  {
    for (j=0; j<NUM_LINEAGES; j++)
    {
      if (species_info->present[j])
      {
        if (species_index == i)
        {
            // We at the required species
            species = j;
        }
        species_index = species_index + 1;
      }
    }    
  }
  map_species_enum_to_str(species, pure_lineage_name);
  return pure_lineage_name->buff;
}


int get_ith_lineage_coverage(SpeciesInfo* species_info,int i)
{
  int covg=0;
  int j;
  int species_index=0;
  if (i > species_info->num_species -1 ){
    die("We only have %i species, we can't find %i \n",  species_info->num_species,i+1);
  }
  else
  {
    for (j=0; j<NUM_LINEAGES; j++)
    {
      if (species_info->present[j])
      {
        if (species_index == i)
        {
            // We at the required species
            covg = species_info->median_coverage[j];
        }
        species_index = species_index + 1;
      }
    }    
  }
  return covg;

}

char* get_pure_lineage_name(SpeciesInfo* species_info)
{
  return get_ith_lineage_name(species_info, 0 );
}

int get_pure_lineage_coverage(SpeciesInfo* species_info)
{
  return get_ith_lineage_coverage(species_info, 0 );
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
      return true;
  }
  else
  {
    return false;
  }
}




void find_which_lineage_panels_are_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels)
{
  int i;
  for (i=0; i<NUM_LINEAGES; i++)
  {
    printf("lineage %d : %d \n", i, percentage_coverage[i]);
    if (is_percentage_coverage_above_threshold(percentage_coverage[i],90))
    {
      *num_panels = *num_panels +1;
    }
  }
}

void find_which_species_panels_are_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels)
{
  int i;
  for (i=0; i<NUM_SPECIES; i++)
  {
    printf("species %d : %d \n", i, percentage_coverage[i]);
    if (is_percentage_coverage_above_threshold(percentage_coverage[i],30))
    {
      *num_panels = *num_panels +1;
    }
  }
}

void are_mtbc_and_ntm_present(int* percentage_coverage,boolean* present, 
                                  int* num_panels)
{
  boolean is_MTBC_present = is_percentage_coverage_above_threshold(percentage_coverage[0],70);
  boolean is_NTM_present = is_percentage_coverage_above_threshold(percentage_coverage[1],70);
  present[0] = is_MTBC_present;
  present[1] = is_NTM_present;
  printf("MTBC : %d \n",percentage_coverage[0] );
  printf("NTM : %d \n",percentage_coverage[1] );


  int i;
  for (i=0; i<2; i++)
  {
    if (present[i])
    {
      *num_panels = *num_panels +1;
    }
  }

}



boolean sample_is_mixed(boolean NTM_is_present,boolean MTBC_is_present)

{
  boolean is_mixed = NTM_is_present && MTBC_is_present;
  return (is_mixed);
}

boolean sample_is_MTBC(boolean NTM_is_present,boolean MTBC_is_present)
{
  boolean is_MTBC = MTBC_is_present && !NTM_is_present;
  return (is_MTBC);
}

boolean sample_is_NTM(boolean NTM_is_present,boolean MTBC_is_present)
{
  boolean is_NTM = !MTBC_is_present && NTM_is_present;
  return (is_NTM);
}

boolean is_NTM_present(boolean* complex_presence, boolean* species_presence)
{
  // Is combined NTM panel present OR any of the NTM species panels
  return (false);
}

boolean is_MTBC_present(boolean* complex_presence, boolean* species_presence)
{
  // Is combined MTBC panel present OR any of the MTBC species panels

  return (false);
}


SampleType get_sample_type(boolean* complex_presence, boolean* species_presence)
{
  boolean NTM_is_present = is_NTM_present(complex_presence,  species_presence);
  boolean MTBC_is_present = is_MTBC_present(complex_presence,  species_presence);


  SampleType sample_type;
  if ( sample_is_mixed(NTM_is_present,MTBC_is_present) )
    {
      sample_type=MixedTB;
    }
  else{
    if ( sample_is_MTBC(NTM_is_present,MTBC_is_present) ){
      sample_type = PureMTBC;
    }
    else if ( sample_is_NTM(NTM_is_present,MTBC_is_present) )
    {
      sample_type = PureNTM;
    }
    else
    {
      sample_type = NonTB;
    }    

  }
  return sample_type;  
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
  panel_file_paths[rhodesiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[rhodesiae], "data/tb/species/rhodesiae.fa");
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


SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{
  // Complex
  int num_complex_panels_present =0 ;
  StrBuf* mtbc_and_ntm_file_paths[NUM_COMPLEX];
  load_all_mtbc_and_ntm_file_paths(mtbc_and_ntm_file_paths,install_dir);
  int mtbc_and_ntm_percentage_coverage[NUM_COMPLEX]; // for storing the percentage coverage of each reference
  int mtbc_and_ntm_median_coverage[NUM_COMPLEX]; //median covg
  boolean mtbc_and_ntm_presence[NUM_COMPLEX];
  get_coverage_on_panels(mtbc_and_ntm_percentage_coverage,mtbc_and_ntm_median_coverage,
                          mtbc_and_ntm_file_paths,max_branch_len,db_graph,
                          ignore_first,ignore_last,NUM_COMPLEX);
  are_mtbc_and_ntm_present(mtbc_and_ntm_percentage_coverage,mtbc_and_ntm_presence,
                                &num_complex_panels_present);
  //species
  int num_species_present =0 ;
  StrBuf* species_file_paths[NUM_SPECIES];
  load_all_species_file_paths(species_file_paths,install_dir);
  int species_percentage_coverage[NUM_SPECIES]; // for storing the percentage coverage of each reference
  int species_median_coverage[NUM_SPECIES]; //median covg
  boolean species_presence[NUM_SPECIES];
  get_coverage_on_panels(species_percentage_coverage,species_median_coverage,
                          species_file_paths,max_branch_len,db_graph,
                          ignore_first,ignore_last,NUM_SPECIES);
  find_which_species_panels_are_present(species_percentage_coverage,species_presence,
                              &num_species_present);
  SampleType sample_type = get_sample_type(mtbc_and_ntm_presence,species_presence);


  // Lineage (ONLY if MTBC - tuberculois)

  int num_lineages_present =0 ;
  StrBuf* lineage_file_paths[NUM_LINEAGES];
  load_all_lineage_file_paths(lineage_file_paths,install_dir);
  int lineage_percentage_coverage[NUM_LINEAGES]; // for storing the percentage coverage of each reference
  int lineage_median_coverage[NUM_LINEAGES]; //median covg
  boolean lineage_presence[NUM_LINEAGES];
  get_coverage_on_panels(lineage_percentage_coverage,lineage_median_coverage,
                          lineage_file_paths,max_branch_len,db_graph,
                          ignore_first,ignore_last,NUM_LINEAGES);
  find_which_lineage_panels_are_present(lineage_percentage_coverage,lineage_presence,
                              &num_lineages_present);


  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
  species_info->sample_type = sample_type;
  species_info->num_species = num_lineages_present;

  memcpy (species_info->present, lineage_presence, sizeof(lineage_presence));
  memcpy (species_info->percentage_coverage, lineage_percentage_coverage, sizeof(lineage_percentage_coverage));
  memcpy (species_info->median_coverage, lineage_median_coverage, sizeof(lineage_median_coverage));
  return species_info;
}




