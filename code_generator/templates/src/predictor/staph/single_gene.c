genotyped_present = false;
I= resistotype_gene(abi->genes[{{drug.genes_resistance_induced_by[0]}}], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon,expected_covg,
         &best_model, MaxAPosteriori,
         {% if drug.name =="Penicillin" %} MIN_PERC_COVG_BLAZ {% else %} MIN_PERC_COVG_STANDARD {% endif %},

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
  if ( genotyped_present ||  cmd_line->verbose) {
    update_called_genes(called_genes, {{drug.genes_resistance_induced_by[0]}}, abi->genes[{{drug.genes_resistance_induced_by[0]}}], best_model.conf );
  }
  update_infection_type(&I,&I_permanent);