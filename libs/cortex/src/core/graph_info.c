/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo  
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  graph_info.c
*/

#include <string.h>
#include <global.h>
#include <stdlib.h>

// cortex_var headers
#include "dB_graph.h"
#include "graph_info.h"
#include "file_reader.h"

int MAX_LEN_SAMPLE_NAME=1000;

ErrorCleaning* error_cleaning_alloc_and_init()
{
  ErrorCleaning* ec = (ErrorCleaning*) malloc(sizeof(ErrorCleaning));
  if (ec==NULL)
    {
      die("Cannot even alloc a tiny error-cleaning info object. Must be some OOM problem. Abort\n");
    }
  memset(ec, 0, sizeof(ErrorCleaning));
  ec->name_of_graph_against_which_was_cleaned = (char*) malloc(sizeof(char)*MAX_FILENAME_LENGTH);
  if (ec->name_of_graph_against_which_was_cleaned==NULL)
    {
      die("Cannot even alloc a tiny error-cleaning info object. Must be some OOM problem. Abort\n");      
    }
  error_cleaning_initialise(ec);
  return ec;
}
void error_cleaning_free(ErrorCleaning* cl)
{
  free(cl->name_of_graph_against_which_was_cleaned);
  free(cl);
}

void error_cleaning_initialise(ErrorCleaning* cl)
{
  error_cleaning_initialise_except_pool_cleaning(cl);

  cl->cleaned_against_another_graph=false;

  set_string_to_null(cl->name_of_graph_against_which_was_cleaned, MAX_FILENAME_LENGTH);
  strcat(cl->name_of_graph_against_which_was_cleaned, "undefined");

  cl->len_name_of_graph_against_which_was_cleaned=
    strlen(cl->name_of_graph_against_which_was_cleaned);
}


void error_cleaning_initialise_except_pool_cleaning(ErrorCleaning* cl)
{
  cl->tip_clipping=false;
  cl->remv_low_cov_sups=false;
  cl->remv_low_cov_nodes=false;
  cl->remv_low_cov_sups_thresh=-1;
  cl->remv_low_cov_nodes_thresh=-1;
}

void error_cleaning_assign_with_OR(ErrorCleaning* target, ErrorCleaning* src, boolean dont_set_pool_cleaning)
{
  target->tip_clipping |= src->tip_clipping;
  target->remv_low_cov_sups |= src->remv_low_cov_sups;
  target->remv_low_cov_nodes |= src->remv_low_cov_nodes;
  target->remv_low_cov_sups_thresh= src->remv_low_cov_sups_thresh;
  target->remv_low_cov_nodes_thresh= src->remv_low_cov_nodes_thresh;
  if ((dont_set_pool_cleaning==false) && (strcmp(src->name_of_graph_against_which_was_cleaned, "") !=0))
    {
      target->cleaned_against_another_graph |=src->cleaned_against_another_graph;
      target->len_name_of_graph_against_which_was_cleaned= src->len_name_of_graph_against_which_was_cleaned;
      target->name_of_graph_against_which_was_cleaned[0]='\0';
      strcat(target->name_of_graph_against_which_was_cleaned, src->name_of_graph_against_which_was_cleaned);
    }

}

GraphInfo* graph_info_alloc_and_init()
{
  GraphInfo* ginfo = (GraphInfo*) malloc(sizeof(GraphInfo));
  if (ginfo==NULL)
    {
      die("Cannot even alloc a GraphInfo object. Your machine must be severely out of memory. Abort\n");
    }
  else
    {
      memset(ginfo, 0, sizeof(GraphInfo));
      int i;
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  ginfo->sample_ids[i] = (char*) malloc(sizeof(char)*MAX_LEN_SAMPLE_NAME) ;
	  if (ginfo->sample_ids[i]==NULL)
	    {
	      die("Cannot even alloc a GraphInfo object. Your machine must be severely out of memory. Abort\n");
	    }
	  ginfo->cleaning[i]=error_cleaning_alloc_and_init(); //will abort if cannot alloc. Should never happen
	}
      
      graph_info_initialise(ginfo);
      return ginfo;
    }
  return NULL;
}

void graph_info_free(GraphInfo* ginfo)
{
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      free(ginfo->sample_ids[i]);
      error_cleaning_free(ginfo->cleaning[i]);
    }
  free(ginfo);
}

void graph_info_initialise(GraphInfo* ginfo)
{
  int i;

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      set_string_to_null(ginfo->sample_ids[i], MAX_LEN_SAMPLE_NAME);
      strcat(ginfo->sample_ids[i], "undefined");
      ginfo->sample_id_lens[i]=strlen(ginfo->sample_ids[i]);
      graph_info_set_seq(ginfo, i, 0);
      graph_info_set_mean_readlen(ginfo, i, 0);
      ginfo->seq_err[i]=0.01;
      error_cleaning_initialise(ginfo->cleaning[i]);
    }
}

void graph_info_initialise_one_colour_except_pool_cleaning(GraphInfo* ginfo, int colour)
{
  set_string_to_null(ginfo->sample_ids[colour], MAX_LEN_SAMPLE_NAME);
  strcat(ginfo->sample_ids[colour], "undefined");
  ginfo->sample_id_lens[colour]=strlen(ginfo->sample_ids[colour]);
  graph_info_set_seq(ginfo, colour, 0);
  graph_info_set_mean_readlen(ginfo, colour, 0);
  ginfo->seq_err[colour]=0.01;
  error_cleaning_initialise_except_pool_cleaning(ginfo->cleaning[colour]);
}





void graph_info_set_tip_clipping(GraphInfo* ginfo, int colour)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->tip_clipping=true;
}

void graph_info_unset_tip_clipping(GraphInfo* ginfo, int colour)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->tip_clipping=false;
}

void graph_info_set_remv_low_cov_sups(GraphInfo* ginfo, int colour, int thresh)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->remv_low_cov_sups=true;
  ginfo->cleaning[colour]->remv_low_cov_sups_thresh = thresh;
}

void graph_info_unset_remv_low_cov_sups(GraphInfo* ginfo, int colour)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->remv_low_cov_sups=false;
  ginfo->cleaning[colour]->remv_low_cov_sups_thresh = -1;
}

void graph_info_set_remv_low_cov_nodes(GraphInfo* ginfo, int colour, int thresh)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->remv_low_cov_nodes=true;
  ginfo->cleaning[colour]->remv_low_cov_nodes_thresh=thresh;
}

void graph_info_unset_remv_low_cov_nodes(GraphInfo* ginfo, int colour)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->remv_low_cov_nodes=false;
  ginfo->cleaning[colour]->remv_low_cov_nodes_thresh=-1;

}

void graph_info_set_seq_err(GraphInfo* ginfo, int colour, long double err)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->seq_err[colour]=err;
}
void graph_info_set_specific_colour_to_cleaned_against_pool(GraphInfo* ginfo, int colour, char* multicol_binary, int colour_in_multicol)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  //set the "name" to be "binary_blah.ctx colour 5", for example
  ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned[0]='\0';
  strcat(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned, multicol_binary);
  strcat(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned, " (colour ");
  char col_as_str[50];
  col_as_str[0]='\0';
  sprintf(col_as_str, "%d", colour_in_multicol);
  strcat(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned, col_as_str);
  strcat(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned, ")");

  //then set the len variable
  ginfo->cleaning[colour]->len_name_of_graph_against_which_was_cleaned
    = strlen(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned);
  
}

void graph_info_unset_specific_colour_from_cleaned_against_pool(GraphInfo* ginfo, int colour)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned[0]='\0';
  strcat(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned, "undefined");
  ginfo->cleaning[colour]->len_name_of_graph_against_which_was_cleaned=
    strlen(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned);
}

void graph_info_set_seq(GraphInfo* ginfo, int colour, long long num_bp)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->total_sequence[colour]=num_bp;
}


long long graph_info_increment_seq(GraphInfo* ginfo, int colour, long long num_bp)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->total_sequence[colour]+=num_bp;
  return ginfo->total_sequence[colour];
}

void graph_info_set_mean_readlen(GraphInfo* ginfo, int colour, int len)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  ginfo->mean_read_length[colour]=len;
}


//note if you are updating both mean read len AND total seq, 
// then do this one first (once you update the total seq, you no longer know 
//  what it used to be, ie the previous_seq)
int graph_info_update_mean_readlen(GraphInfo* ginfo, int colour, int previous_mean, long long previous_seq, 
			int mean_readlen_in_added_data, long long added_seq)
{
  if (colour>NUMBER_OF_COLOURS-1)
    {
      die("Setting a graph_info error cleaning flag with a colour (%d) that is bigger than you have compiled for (%d). Coding error, the UI should have prevented this. Call Zam\n", colour, NUMBER_OF_COLOURS);
    }
  if (added_seq==0)
    {
      return previous_mean;
    }
  else if (mean_readlen_in_added_data==0)
    {
      //printf("Warning - adding data with mean read-length 0\n");
      return previous_mean;
    }
  long long numerator = ((long long) previous_mean) * previous_seq + 
                         ((long long)mean_readlen_in_added_data) * added_seq;

  if (previous_seq+added_seq==0)
    {
      printf("WARNING - Updating graph_info object which contains no data with extra zero data!\n");
      return 0;
    }
  int new_mean = (int) (numerator/( previous_seq+added_seq));
  ginfo->mean_read_length[colour]=new_mean;
  return new_mean;
}

void graph_info_update_mean_readlen_and_total_seq(GraphInfo* ginfo, int colour,
                                                  unsigned long mean_readlen_in_added_data,
						  unsigned long long added_seq)
{
  graph_info_update_mean_readlen(ginfo, colour, 
		      ginfo->mean_read_length[colour], ginfo->total_sequence[colour], 
		      mean_readlen_in_added_data, added_seq);
  graph_info_increment_seq(ginfo, colour, added_seq);
}


void graph_info_set_sample_ids(char** sample_ids, int num_ids, GraphInfo* ginfo, int first_col)
{
  if (first_col+num_ids-1>NUMBER_OF_COLOURS-1)
    {
      die("Coding error in graph_info_set_sample_ids - trying to load sample id into a colour greater than the number of supported colours you have compiled for. Should have been caught earlier - please tell Zam\n");
    }
  int i;
  for (i=first_col; i<first_col+num_ids; i++)
    {
      strcpy(ginfo->sample_ids[i], sample_ids[i]);
      ginfo->sample_id_lens[i]=strlen(ginfo->sample_ids[i]);
    }
}

//this is designed to be used when loading multiple binaries into a 
//single colour, which may have different metadata. In this case,
// we cumulate covg, get mean read length, and for each type of error cleaning we 
//say yesif ANY of them did, taking the last threshold for any of them.
void graph_info_set_all_metadata(GraphInfo* target, GraphInfo* src, int colour, boolean dont_set_pool_cleaning)
{
  graph_info_update_mean_readlen(target, colour, 
				 target->mean_read_length[colour], 
				 target->total_sequence[colour],
				 src->mean_read_length[colour],
				 src->total_sequence[colour]);
  graph_info_increment_seq(target, colour, src->total_sequence[colour]);
  target->seq_err[colour]=src->seq_err[colour];
  error_cleaning_assign_with_OR(target->cleaning[colour], src->cleaning[colour], dont_set_pool_cleaning);


  if ( (strcmp(src->sample_ids[colour], "undefined")!=0) //source has nontrivial sample id
       &&
       (strcmp(src->sample_ids[colour], target->sample_ids[colour])!=0) //src and target have different sample id
       && 
       (strcmp(target->sample_ids[colour],"undefined")!=0)//target has non trivial sample id
       )
    {
      set_string_to_null(target->sample_ids[colour], MAX_LEN_SAMPLE_NAME );
      strcat(target->sample_ids[colour], "pool");
    }
  else if (strcmp(src->sample_ids[colour], "undefined")!=0) //src has nontrivial sample id
    {
      //either it is the same, or target is undefined. Either way, this is ok:
      set_string_to_null(target->sample_ids[colour], MAX_LEN_SAMPLE_NAME);
      strcat(target->sample_ids[colour], src->sample_ids[colour]);
    }
}

double get_total_depth_of_coverage_across_colours(GraphInfo* ginfo, long long genome_length)
{
  int i;
  long long total_seq =0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      total_seq += ginfo->total_sequence[i];
    }
  
  if (total_seq/genome_length<1)
    {
      printf("Warning - total sequence contained in union of all colours is less than 1x coverage of your genome. Is this really what you intend?\n");
    }

  return ((double) total_seq)/ ((double) genome_length);
}


int get_mean_readlen_across_colours(GraphInfo* ginfo)
{
  int colour;

  long long alpha=0;
  long long beta=0;
  for (colour=0; colour<NUMBER_OF_COLOURS; colour++)
    {
      alpha += (ginfo->total_sequence[colour]) * (ginfo->mean_read_length[colour]);
      beta += (ginfo->total_sequence[colour]);
    }
  return (int) (alpha/beta);
}

void read_estimated_seq_errors_from_file(GraphInfo* ginfo, FILE* fp)
{
  char line[MAX_FILENAME_LENGTH+1];
  int col=0;

  while (fgets(line, MAX_FILENAME_LENGTH, fp) !=NULL)
    {
            //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      if (col<NUMBER_OF_COLOURS) //just for robustness - should always be true
	{
	  ginfo->seq_err[col]=  (long double) strtod(line ,NULL);
	  col++;
	}
	
    }
}

void print_seq_err_rates_to_screen(GraphInfo* ginfo)
{
  printf("Setting the following per-colour sequencing error rates (used only for genotyping):\nColour\tRate\n");
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      printf("%d\t%.3Lf\n", i, ginfo->seq_err[i]);
    }
}
