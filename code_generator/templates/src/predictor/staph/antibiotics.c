{% extends 'src/predictor/common/antibiotics.c' %}

{% block additional_functions %}

void print_clindamycin_susceptibility(dBGraph* db_graph,
					 int (*file_reader)(FILE * fp, 
							    Sequence * seq, 
							    int max_read_length, 
							    boolean new_entry, 
							    boolean * full_entry),
					 ReadingUtils* rutils,
					 VarOnBackground* tmp_vob,
					 GeneInfo* tmp_gi,
					 AntibioticInfo* abi,
					 InfectionType (*func)(dBGraph* db_graph,
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
               CalledVariant* called_variants,CalledGene* called_genes,
               CmdLine* cmd_line),
					 StrBuf* tmpbuf,
					 boolean any_erm_present, InfectionType erythromycin_resistotype,
					 StrBuf* install_dir,
					 int ignore_first, int ignore_last, SpeciesInfo* species_info,
					 double lambda_g, double lambda_e, double err_rate, CmdLine* cmd_line, boolean output_last,
           CalledVariant* called_variants,CalledGene* called_genes//for JSON 
					 )
{
  InfectionType suc;
  
  suc  = func(db_graph,
	      file_reader,
	      rutils,
	      tmp_vob,
	      tmp_gi,
	      abi,
	      install_dir,
	      ignore_first, ignore_last, species_info,
	      lambda_g,
	      lambda_e,
	      err_rate,
        called_variants,
        called_genes,
        cmd_line);


  map_antibiotic_enum_to_str(abi->ab, tmpbuf);


    if (suc==Resistant)
	{
	  print_json_item(tmpbuf->buff, "R(constitutive)", output_last);
	}
      else if (suc==MixedInfection)
	{
	  print_json_item(tmpbuf->buff, "r(constitutive)", output_last);
	}
      else if ( (suc==Susceptible) && (any_erm_present==true) )
	{
	    if (erythromycin_resistotype == Resistant){
	      print_json_item(tmpbuf->buff, "R(inducible)", output_last);
	    }
	    else if (erythromycin_resistotype == MixedInfection){
	      print_json_item(tmpbuf->buff, "r(inducible)", output_last);
	    }	
	}
      else if (suc==Susceptible)
	{
	  print_json_item(tmpbuf->buff, "S", output_last);
	}
      else
	{
	  print_json_item(tmpbuf->buff, "Inconclusive", output_last);
	}

}


void print_erythromycin_susceptibility(dBGraph* db_graph,
					  int (*file_reader)(FILE * fp, 
							     Sequence * seq, 
							     int max_read_length, 
							     boolean new_entry, 
							     boolean * full_entry),
					  ReadingUtils* rutils,
					  VarOnBackground* tmp_vob,
					  GeneInfo* tmp_gi,
					  AntibioticInfo* abi,
					  InfectionType (*func)(dBGraph* db_graph,
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
								boolean* any_erm_present, 
                CalledVariant* called_variants,CalledGene* called_genes,
                CmdLine* cmd_line),
					  StrBuf* tmpbuf,
					  StrBuf* install_dir,
					  int ignore_first, int ignore_last, SpeciesInfo* species_info,
					  double lambda_g, double lambda_e, double err_rate, CmdLine* cmd_line, boolean output_last,//for JSON 
					  boolean* any_erm_present, InfectionType* erythromycin_resistotype, 
            CalledVariant* called_variants,CalledGene* called_genes
					 )
{
  InfectionType suc;
  
  suc  = func(db_graph,
	      file_reader,
	      rutils,
	      tmp_vob,
	      tmp_gi,
	      abi,
	      install_dir,
	      ignore_first, ignore_last, species_info,
	      lambda_g,
	      lambda_e,
	      err_rate,
	      any_erm_present,
        called_variants,
         called_genes,
         cmd_line);

  map_antibiotic_enum_to_str(abi->ab, tmpbuf);
  *erythromycin_resistotype = suc;
  
  if (cmd_line->format==Stdout)
    {
      printf("%s\t", tmpbuf->buff);
      if (suc==Susceptible)
	{
	  printf("S\n");
	}
      else if (suc==MixedInfection)
	{
	  printf("r\n");
	}
      else if (suc==Resistant)
	{
	  printf("R\n");
	}
      else
	{
	  printf("N\n");
	}
    }
  else
    {
      if (suc==Susceptible)
	{
	  print_json_item(tmpbuf->buff, "S", output_last);
	}
      else if ( suc==Resistant)
	{
	  print_json_item(tmpbuf->buff, "R", output_last);
	}
      else if ( suc==MixedInfection ) 
  {
    print_json_item(tmpbuf->buff, "r", output_last);
  }
      else
	{
	  print_json_item(tmpbuf->buff, "Inconclusive", output_last);
	}
    }

}


{% endblock %}
