  int mut = {{drug.first_mut}};
  i = mut;
  boolean any_allele_non_null=false;
  if (both_alleles_null(abi->vars[i])==true)
	{
	  any_allele_non_null=true;
	}
  genotyped_present = false;
  InfectionType I=resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
    lambda_g, lambda_e, epsilon, expected_covg, 
    &best_model, MaxAPosteriori,
    cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);
  if ( genotyped_present ||  cmd_line->verbose) {
    update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  } 

  if ( (I==Susceptible) && (best_model.conf>max_sus_conf) )
	{
	  max_sus_conf = best_model.conf;
	}
  if (best_model.conf<min_conf)
	{
	  min_conf = best_model.conf;
	}
  update_infection_type(&I,&I_permanent);