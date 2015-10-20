InfectionType is_{{drug | lower }}_susceptible(dBGraph* db_graph,
				  int (*file_reader)(FILE * fp, 
						     Sequence * seq, 
						     int max_read_length, 
						     boolean new_entry, 
						     boolean * full_entry),
				  ReadingUtils* rutils,
				  VarOnBackground* tmp_vob,
				  GeneInfo* tmp_gi,
				  AntibioticInfo* abi,
				  StrBuf* install_dir,
				  int ignore_first, int ignore_last, SpeciesInfo* species_info,
				  double lambda_g, double lambda_e, double err_rate,
            	  {% if drug.name =="Erythromycin" %} boolean* any_erm_present,{% endif %}				  
				   CalledVariant* called_variants,CalledGene* called_genes,
				   CmdLine* cmd_line
				  );