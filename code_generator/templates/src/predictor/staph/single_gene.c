  InfectionType I=
    resistotype_gene(abi->genes[{{drug.genes_resistance_induced_by[0]}}], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon,expected_covg,
         &best_model, MaxAPosteriori,
         MIN_PERC_COVG_BLAZ);
  if ( (I==Resistant) || (I==MixedInfection) ) {
    update_called_genes(called_genes, {{drug.genes_resistance_induced_by[0]}}, abi->genes[{{drug.genes_resistance_induced_by[0]}}], best_model.conf );
  }
  update_infection_type(&I,&I_permenant);