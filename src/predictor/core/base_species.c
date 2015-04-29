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
#include "base_species.h"
#include "gene_presence.h"
#include "genotyping_known.h"


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
    covg_info->present[i] = is_percentage_coverage_above_threshold(covg_info->percentage_coverage[i],covg_info->percentage_threshold[i]);
    if (covg_info->present[i])
    {
      covg_info->num_panels_present = covg_info->num_panels_present +1;
    }
  }
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

int max(int int1,int int2){
  if (int1 > int2){
    return (int1);
  }
  else{
    return (int2);
  }

}

int get_best_hit(CovgInfo* covg_info,boolean* mask)
{
  int i;
  int best_perc_cov_so_far=0;
  int best_median_cov_so_far=0;
  Species unknown_enum = unknown;
  int curr=(int) unknown_enum;

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
  return  curr;
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

char* get_char_name_of_species_enum(Species species){
  StrBuf* species_name = strbuf_new(); 
  map_species_enum_to_str(species, species_name);
  return species_name->buff;
}



char* get_ith_phylo_group_name(CovgInfo* covg_info, int i)
{
  PhyloGroup phylo_group;
  StrBuf* phylo_group_name = strbuf_new(); 
  phylo_group = get_ith_present_panel( covg_info, i);
  map_phylo_group_enum_to_str(phylo_group, phylo_group_name);
  return phylo_group_name->buff;
}
char* get_ith_species_name(CovgInfo* covg_info, int i)
{
  Species species;
  species = get_ith_present_panel( covg_info, i);
  return get_char_name_of_species_enum(species);
}



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

void print_json_phylogenetics(SpeciesInfo* species_info){
    print_json_phylogenetics_start();
    print_json_phylo_group(species_info);
    print_json_species(species_info);
    print_json_lineage(species_info);
    print_json_phylogenetics_end();  
}

#ifdef TB

char* get_char_name_of_lineage_enum(Lineage lineage){
  StrBuf* lineage_name = strbuf_new(); 
  map_lineage_enum_to_str(lineage, lineage_name);
  return lineage_name->buff;
}

char* get_ith_lineage_name(CovgInfo* covg_info, int i)
{
  Lineage lineage;
  lineage = get_ith_present_panel( covg_info, i);
  return get_char_name_of_lineage_enum(lineage);
}

 #endif
