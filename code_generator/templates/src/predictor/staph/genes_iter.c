
for (i=0; i<{{drug.num_genes}}; i++)
    {
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg,
			 &best_model, MaxAPosteriori,
			 MIN_PERC_COVG_STANDARD);
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
        {% if drug.name =="Erythromycin" %} 
          if ( (abi->which_genes[i] == ermA) || (abi->which_genes[i] == ermB) || (abi->which_genes[i] == ermC) || (abi->which_genes[i] == ermY) || (abi->which_genes[i] == ermT) )
          {
            *any_erm_present=true;
          }
        {% endif %}
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
         else if (cmd_line->verbose){
           update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
         }        
        update_infection_type(&I,&I_permenant);
    }