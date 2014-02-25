
// cortex_var headers
#include "db_variants.h"
#include "maths.h"


char variant_overflow_warning_printed = 0;




//first argument is an array of length NUMBER_OF_COLOURS, into which results go.
//If you want the number of reads on the entire branch, enter the length of that branch in arg3 (eg var->len_one-allele)
//Sometimes we want to take just the start of the branch (if one branch is longer than the other, we may just take the length of the shorter one)
//and so you enter that in arg3 in that case
//note these are effective reads, as counting covg in the de Bruijn graph
//returns TRUE if branch is too short (1 oor 2 nodes) to do this
boolean get_num_effective_reads_on_branch(Covg* array, dBNode** allele, int how_many_nodes, 
					  boolean use_median, CovgArray* working_ca, GraphInfo* ginfo, int kmer)
{
  int i;
  boolean too_short=false;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (use_median==false)
	{
	  array[i] = count_reads_on_allele_in_specific_colour(allele, how_many_nodes, i, &too_short);
	}
      else
	{
	  int eff_read_len = 100;
	  if (ginfo!=NULL)
	    {
	      eff_read_len = ginfo->mean_read_length[i] - kmer+1;
	    }
	  if (how_many_nodes>eff_read_len)
	    {
	      array[i] = ((how_many_nodes+ 0.5*eff_read_len)/eff_read_len) *  median_covg_on_allele_in_specific_colour(allele, how_many_nodes, working_ca, i, &too_short);
	    }
	  else
	    {
	      array[i] = median_covg_on_allele_in_specific_colour(allele, how_many_nodes, working_ca, i, &too_short);
	    }
	}
    }
  return too_short;
}






//does not count covg on first or last nodes, as they are bifurcation nodes
//if length==0 or 1  returns 0.
//note I do not want to create an array on the stack - these things can be very long
// so relies on the prealloced array of dBNode* 's passed in
// annoying that can't use templates or something - see below for a similar function with different input
Covg count_reads_on_allele_in_specific_colour(dBNode** allele, int len, int colour,
                                              boolean* too_short)
{
  if(len <= 1)
  {
    *too_short = true;
    return 0;
  }

  // Note: we have to used signed datatypes for this arithmetic
  //       hence using (long) instead of (Covg -- usually uint32_t)

  // note start at node 1, avoid first node
  Covg num_of_reads = db_node_get_coverage_tolerate_null(allele[1], colour);

  int i;

  //note we do not go as far as the final node, which is where the two branches rejoin
  for(i = 2; i < len-1; i++)
  {
    long jump = (long)db_node_get_coverage_tolerate_null(allele[i], colour) -
                (long)db_node_get_coverage_tolerate_null(allele[i-1], colour);

    // we add a little check to ensure that we ignore isolated nodes with higher
    // covg - we count jumps only if they are signifiers of a new read arriving
    // and one node does not a read make
    long diff_between_next_and_prev = -1;

    if(i < len-2)
  	{
      diff_between_next_and_prev
        = (long)db_node_get_coverage_tolerate_null(allele[i+1], colour) -
          (long)db_node_get_coverage_tolerate_null(allele[i-1], colour);
    }

    if(jump > 0 && diff_between_next_and_prev != 0)
    {
      if(COVG_MAX - jump >= num_of_reads)
      {
        num_of_reads += jump;
      }
      else
      {
        num_of_reads = COVG_MAX;

        if(!variant_overflow_warning_printed)
        {
          warn("%s:%i: caught integer overflow"
               "(some kmer coverages may be underestimates)",
               __FILE__, __LINE__);

          variant_overflow_warning_printed = 1;
        }
      }
    }
  }

  return num_of_reads;
}



// WARNING - this is for use when we dissect an allele into subchunks, so here
// we do not want to be ignoring first and last elements (cf above)
Covg count_reads_on_allele_in_specific_colour_given_array_of_cvgs(Covg* covgs,
                                                                  int len,
                                                                  boolean* too_short)
{
  if(len <= 1)
  {
    *too_short = true;
    return 0;
  }

  // Note: we have to used signed datatypes for this arithmetic
  //       hence using (long) instead of (Covg -- usually uint32_t)

  Covg num_of_reads = covgs[0];

  int i;

  for(i = 1; i < len; i++)
  {
    long jump = (long)covgs[i] - covgs[i-1];

    // we add a little check to ensure that we ignore isolated nodes with higher
    // covg - we count jumps only if they are signifiers of a new read arriving
    // and one node does not a read make
    long diff_between_next_and_prev = -1;

    if(i < len-1)
    {
      diff_between_next_and_prev = (long)covgs[i+1] - covgs[i-1];
    }

    if(jump > 0 && diff_between_next_and_prev != 0)
    {
      if(COVG_MAX - jump >= num_of_reads)
      {
        num_of_reads += jump;
      }
      else
      {
        num_of_reads = COVG_MAX;

        if(!variant_overflow_warning_printed)
        {
          warn("%s:%i: caught integer overflow"
               "(some kmer coverages may be underestimates)",
               __FILE__, __LINE__);

          variant_overflow_warning_printed = 1;
        }
      }
    }
  }

  return num_of_reads;
}


//does not count covg on first or last nodes, as they are bifurcation nodes
//if length==0 or 1  returns 0.
Covg count_reads_on_allele_in_specific_func_of_colours(
  dBNode** allele, int len,
  Covg (*sum_of_covgs_in_desired_colours)(const Element *),
  boolean* too_short)
{
  if(len <= 1)
  {
    *too_short = true;
    return 0;
  }

  // Note: we have to used signed datatypes for this arithmetic
  //       hence using (long) instead of (Covg -- usually uint32_t)

  // note start at node 1, avoid first node
  Covg num_of_reads = sum_of_covgs_in_desired_colours(allele[1]);

  int i;

  //note we do not go as far as the final node, which is where the two branches rejoin
  for(i = 2; i < len-1; i++)
  {
    long jump = (long)sum_of_covgs_in_desired_colours(allele[i]) -
                (long)sum_of_covgs_in_desired_colours(allele[i-1]);

    // we add a little check to ensure that we ignore isolated nodes with higher
    // covg - we count jumps only if they are signifiers of a new read arriving
    // and one node does not a read make
    long diff_between_next_and_prev = -1;
    
    if(i < len-2)
    {
      diff_between_next_and_prev
        = (long)sum_of_covgs_in_desired_colours(allele[i+1]) -
          (long)sum_of_covgs_in_desired_colours(allele[i-1]);
    }
      
    if(jump > 0 && diff_between_next_and_prev != 0)
    {
      // Increment number of reads -- avoid overflow
      if(COVG_MAX - jump >= num_of_reads)
      {
        num_of_reads += jump;
      }
      else
      {
        num_of_reads = COVG_MAX;

        if(!variant_overflow_warning_printed)
        {
          warn("%s:%i: caught integer overflow"
               "(some kmer coverages may be underestimates)",
               __FILE__, __LINE__);

          variant_overflow_warning_printed = 1;
        }
      }
    }
  }

  return num_of_reads;
}




//robust to start being > end (might traverse an allele backwards)
//if length==0 or 1  returns 0.
Covg median_covg_on_allele_in_specific_colour(dBNode** allele, int len, CovgArray* working_ca,
					      int colour, boolean* too_short)
{
  printf("Len is %d\n", len);

  if ((len==0)|| (len==1))
    {
      *too_short=true;
      return 0;//ignore first and last nodes
    }
 
  reset_covg_array(working_ca);//TODO - use reset_used_part_of.... as performance improvement. Will do when correctness of method established.
  int i;

  for(i=1; i <len; i++)
    {
      working_ca->covgs[i-1]=db_node_get_coverage_tolerate_null(allele[i], colour);
    }

  int array_len =len-1;
  qsort(working_ca->covgs, array_len, sizeof(Covg), Covg_cmp); 
  working_ca->len=len-1;

  Covg median=0;
  int lhs = (array_len - 1) / 2 ;
  int rhs = array_len / 2 ;
  
  if (lhs == rhs)
    {
      median = working_ca->covgs[lhs] ;
    }
  else 
    {
      median = mean_of_covgs(working_ca->covgs[lhs], working_ca->covgs[rhs]);
    }

  return median;
}



//robust to start being > end (might traverse an allele backwards)
//if length==0 or 1  returns 0.
Covg min_covg_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short)
{

  if ((len==0)|| (len==1))
    {
      *too_short=true;
      return 0;//ignore first and last nodes
    }
 
  int i;

  int index=0;

  Covg min_covg = COVG_MAX;
  for(i=1; i <len; i++)
    {
      if (allele[i]!=NULL)
	{
	  Covg c=db_node_get_coverage(allele[i], colour);
	  if (c<min_covg)
	    {
	      min_covg = c;
	    }
	}
    }

  return min_covg;
}


int percent_nonzero_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short)
{

  if ((len==0)|| (len==1))
    {
      *too_short=true;
      return 0;//ignore first and last nodes
    }
 
  int i;

  int index=0;
  int num_nonzero=0;

  for(i=1; i <len; i++)
    {
      if (allele[i]!=NULL)
	{
	  Covg c=db_node_get_coverage(allele[i], colour);
	  if (c>0)
	    {
	      num_nonzero++;
	    }
	}
    }

  return num_nonzero*100/(len-1);
}


//only count nodes which have the desired allele status
Covg median_covg_on_allele_in_specific_colour_with_allele_presence_constraint(dBNode** allele, int len, CovgArray* working_ca,
									      int colour, boolean* too_short, AlleleStatus st,
									      float eff_depth)
{

  if ((len==0)|| (len==1))
    {
      *too_short=true;
      return 0;//ignore first and last nodes
    }
 
  reset_covg_array(working_ca);//TODO - use reset_used_part_of.... as performance improvement. Will do when correctness of method established.
  int i;

  int index=0;

  for(i=1; i <len; i++)
    {
      if (allele[i]!=NULL)
	{
	  if (db_node_check_allele_status(allele[i], st)==true)
	    {
	      Covg cov = db_node_get_coverage_tolerate_null(allele[i], colour);
	      if (cov < 2* eff_depth)
		{
		  working_ca->covgs[index]=cov;
		  index++;
		}
	    }
	}
    }

  int array_len = index+1;
  qsort(working_ca->covgs, array_len, sizeof(Covg), Covg_cmp); 
  working_ca->len=index+1;

  Covg median=0;
  int lhs = (array_len - 1) / 2 ;
  int rhs = array_len / 2 ;
  
  if (lhs == rhs)
    {
      median = working_ca->covgs[lhs] ;
    }
  else 
    {
      median = mean_of_covgs(working_ca->covgs[lhs], working_ca->covgs[rhs]);
    }

  return median;
}





Covg median_of_CovgArray(CovgArray* array, CovgArray* working_array)
{
  if (array->len>working_array->len_alloced)
    {
      die("Trying to find median of an array using a working-array which is too short\n");
    }
  int i=0;
  for(i = 0; i<array->len; i++)
    {
      working_array->covgs[i]=array->covgs[i];
    }

  qsort(working_array->covgs, array->len, sizeof(Covg), Covg_cmp); 
  
  Covg median=0;
  int lhs = (array->len - 1) / 2 ;
  int rhs = array->len / 2 ;
  
  if (lhs == rhs)
    {
      median = working_array->covgs[lhs] ;
    }
  else 
    {
      median = mean_of_covgs(working_array->covgs[lhs], working_array->covgs[rhs]);
    }
  return median;
}

