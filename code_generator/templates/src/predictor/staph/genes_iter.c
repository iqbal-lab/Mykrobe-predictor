
for (i=0; i<{{drug.num_genes}}; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg,
			 &best_model, MaxAPosteriori,
			 {% if drug.name =="Trimethoprim" %} MIN_PERC_COVG_DFRK {% else %} MIN_PERC_COVG_STANDARD {% endif %},
    
    {% if drug.name =="Penicillin" %} MIN_GENE_CN_PEN
    {% elif drug.name == "Gentamicin" %} MIN_GENE_CN_GEN 
    {% elif drug.name == "Methicillin" %} MIN_GENE_CN_MEC 
    {% elif drug.name == "FusidicAcid" %} MIN_GENE_CN_FUS 
    {% elif drug.name == "Erythromycin" %} MIN_GENE_CN_ERY 
    {% elif drug.name == "Mupirocin" %} MIN_GENE_CN_MUP 
    {% elif drug.name == "Tetracycline" %} MIN_GENE_CN_TET 
    {% elif drug.name == "Trimethoprim" %} MIN_GENE_CN_PEN 
    {% else %} MIN_GENE_CN 
    {% endif %},&genotyped_present

        );
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
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permenant);
    }