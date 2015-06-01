{% if drug.name == "FusidicAcid" %}
  InfectionType I_f652s=
  resistotype(abi->vars[fusA_F652S],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops);

  InfectionType I_y654n=
    resistotype(abi->vars[fusA_Y654N],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops);
if (I_f652s==Resistant && I_y654n==Resistant)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_F652S], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_Y654N], best_model.conf);
    update_infection_type(Resistant,&I_permenant);  
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_F652S], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_Y654N], best_model.conf);    
  }




  InfectionType I_t326i=
  resistotype(abi->vars[fusA_T326I],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops);

  InfectionType I_e468v=
    resistotype(abi->vars[fusA_E468V],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops);


if (I_t326i==Resistant && I_e468v==Resistant)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_T326I],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E468V],best_model.conf);
    update_infection_type(Resistant,&I_permenant);  
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_T326I],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E468V],best_model.conf);
  }

  



  InfectionType I_l461f=
  resistotype(abi->vars[fusA_L461F],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops);

  InfectionType I_a376v=
    resistotype(abi->vars[fusA_A376V],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops);

  InfectionType I_a655p=
    resistotype(abi->vars[fusA_A655P],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops);

  InfectionType I_d463g=
    resistotype(abi->vars[fusA_D463G],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops);
  
if ( (I_l461f==Resistant)
       &&
     (I_a376v==Resistant)
       &&
     (I_a655p==Resistant)
       &&
     (I_d463g==Resistant)
     )
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_A376V],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_A655P],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_D463G],best_model.conf);            
    update_infection_type(Resistant,&I_permenant);  
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_A376V],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_A655P],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_D463G],best_model.conf);
  }

InfectionType I_e444v=Susceptible;
resistotype(abi->vars[fusA_E444V],
	    err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
	    &best_model, MaxAPosteriori,
	    cmd_line->min_frac_to_detect_minor_pops);
if ((I_l461f==Resistant)
       &&
    (I_e444v==Resistant) )
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E444V],best_model.conf); 
    update_infection_type(Resistant,&I_permenant);     
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E444V],best_model.conf); 
  }  
{% elif drug.name =="Rifampicin" %}
  InfectionType I_m470t=
  resistotype(abi->vars[rpoB_M470T],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops);

if ( (I_m470t==Susceptible) && (best_model.conf>max_sus_conf) )
  {
    max_sus_conf = best_model.conf;
  }
if (best_model.conf<min_conf)
  {
    min_conf = best_model.conf;
  }

  InfectionType I_d471g=
  resistotype(abi->vars[rpoB_D471G],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops);

if ( (I_d471g ==Susceptible) && (best_model.conf>max_sus_conf) )
  {
    max_sus_conf = best_model.conf;
  }
if (best_model.conf<min_conf)
  {
    min_conf = best_model.conf;
  }

if (I_m470t==Resistant && I_d471g==Resistant)
  {
    update_called_variants(called_variants,i,abi->vars[rpoB_D471G], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[rpoB_M470T], best_model.conf);
    update_infection_type(Resistant,&I_permenant);   //ignoring mixed infections for epistatic case
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[rpoB_D471G], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[rpoB_M470T], best_model.conf);
  }    
  
{% endif %}
