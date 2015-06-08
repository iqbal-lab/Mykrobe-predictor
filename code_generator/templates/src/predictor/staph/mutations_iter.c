  int first_mut = {{drug.first_mut}};
  int last_mut = {{drug.last_mut}};
  //we are going to iterate through various mutations, each
  //on different genetic backgrounds.
  //we will call it susceptible, if the best hit is good enough
  //if you have any of these resistance alleles - call resistant
  boolean any_allele_non_null=false;
  for (i=first_mut; i<=last_mut; i++)
    {
      if (both_alleles_null(abi->vars[i])==true)
    	{
    	  continue;
    	}
      any_allele_non_null=true;
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,expected_depth, 
		    &best_model, MaxAPosteriori,
		    cmd_line->min_frac_to_detect_minor_pops);
      if ( (I==Susceptible) && (best_model.conf>max_sus_conf) )
    	{
    	  max_sus_conf = best_model.conf;
    	}
      if (best_model.conf<min_conf)
    	{
    	  min_conf = best_model.conf;
    	}
      if ( (I==Resistant) || (I==MixedInfection) ) 
    	{
    	  update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
    	}
      else if (cmd_line->verbose)
      {
        update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
      }  
      update_infection_type(&I,&I_permenant);
    }