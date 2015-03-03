/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
  antibiotics.c 
*/

#include "antibiotics.h"
#include "file_reader.h"
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include "mut_models.h"
#include "gene_presence_models.h"
#include "json.h"

char* map_gene_to_drug_resistance(GenePresenceGene gene)
{
   return "unknown";
}

void map_antibiotic_enum_to_str(Antibiotic ab, StrBuf* name)
{
        if (ab==NoDrug)
        {
          strbuf_reset(name);
          strbuf_append_str(name, "NoDrug");
        }
        else if (ab==rifampicin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Rifampicin");
                }
        
  else if (ab==isoniazid)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Isoniazid");
                }
        
  else if (ab==pyrazinamide)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Pyrazinamide");
                }
        
  else if (ab==ethambutol)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Ethambutol");
                }
        
  else if (ab==kanamycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Kanamycin");
                }
        
  else if (ab==capreomycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Capreomycin");
                }
        
  else if (ab==amikacin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Amikacin");
                }
        
  else if (ab==kanamycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Kanamycin");
                }
        
  else if (ab==streptomycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Streptomycin");
                }
        
  else if (ab==quinolones)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "Quinolones");
                }
         
  else
  {
    die("Impossible - compiler should not allow this");
  }
}

AntibioticInfo* alloc_antibiotic_info()
{
  AntibioticInfo* abi = (AntibioticInfo*) calloc(1, sizeof(AntibioticInfo));
  if (abi==NULL)
    {
      return NULL;
    }
  else
    {
      abi->m_fasta = strbuf_new();
      abi->num_genes = 0;
      abi->vars = (Var**) malloc(sizeof(Var*)*NUM_KNOWN_MUTATIONS);
      if (abi->vars==NULL)
	{
	  strbuf_free(abi->m_fasta);
          free(abi);
	  return NULL;
	}
      abi->genes = (GeneInfo**) malloc(sizeof(GeneInfo*)*NUM_GENE_PRESENCE_GENES);
      if (abi->genes==NULL)
	{
	  free(abi->vars);
	  strbuf_free(abi->m_fasta);
	  free(abi);
	  return NULL;
	}
      abi->which_genes = (int*) calloc(MAX_GENES_PER_ANTIBIOTIC, sizeof(int));
      if (abi->which_genes==NULL)
	{
	  free(abi->genes);
	  free(abi->vars);
	  strbuf_free(abi->m_fasta);
	  free(abi);
	  return NULL;
	  
	}
      int i;
      for (i=0; i<NUM_KNOWN_MUTATIONS; i++)
	{
	  abi->vars[i] = alloc_var();
	  if (abi->vars[i]==NULL)
	    {
	      free(abi->vars);
	      free(abi->genes); 
	      strbuf_free(abi->m_fasta);
	      free(abi);
	      return NULL; //creates a leak if i>0
	    }
	}
      for (i=0; i<NUM_GENE_PRESENCE_GENES; i++)
	{
	  abi->genes[i] = alloc_and_init_gene_info();
	  if (abi->genes[i]==NULL)
	    {
	      free(abi->vars);
              free(abi->genes);
              strbuf_free(abi->m_fasta);
              free(abi);
              return NULL; 
	    }
	  abi->genes[i]->name = (GenePresenceGene) i;
	}

    }
  return abi;
}
void free_antibiotic_info(AntibioticInfo* abi)
{
  if (abi==NULL)
    {
      return;
    }
  else
    {
      strbuf_free(abi->m_fasta);
      int i;
      for (i=0; i<NUM_KNOWN_MUTATIONS; i++)
	{
	  free_var(abi->vars[i]);
	}
      free(abi->vars);
      for (i=0; i<NUM_GENE_PRESENCE_GENES; i++)
	{
	  free_gene_info(abi->genes[i]);
	}
      free(abi->genes);
      free(abi->which_genes);
      free(abi);
    }
}

void reset_antibiotic_info(AntibioticInfo* abi)
{
  abi->ab = NoDrug;
  strbuf_reset(abi->m_fasta);
  abi->num_mutations=0;

  int i;
  for (i=0; i<NUM_KNOWN_MUTATIONS; i++)
    {
      reset_var_on_background(abi->vars[i]->vob_best_sus);
      reset_var_on_background(abi->vars[i]->vob_best_res);
    }
  for (i=0; i<NUM_GENE_PRESENCE_GENES; i++)
    {
      reset_gene_info(abi->genes[i]);
    }
}

void  load_antibiotic_mutation_info_on_sample(FILE* fp,
					      dBGraph* db_graph,
					      int (*file_reader)(FILE * fp, 
								 Sequence * seq, 
								 int max_read_length, 
								 boolean new_entry, 
								 boolean * full_entry),
					      AntibioticInfo* abi,
					      ReadingUtils* rutils,
					      VarOnBackground* tmp_vob,	
					      int ignore_first, int ignore_last, int expected_covg)
{
  reset_reading_utils(rutils);
  reset_var_on_background(tmp_vob);

  StrBuf* tmp1 = strbuf_new();
  StrBuf* tmp2 = strbuf_new();
  StrBuf* tmp3 = strbuf_new();

  
  KnownMutation m = NotSpecified;
  boolean ret=true;

  while (ret==true)
    {
      ret = get_next_var_on_background(fp, 
				       db_graph, 
				       tmp_vob, abi->vars,
				       rutils->seq, 
				       rutils->kmer_window, 
				       file_reader,
				       rutils->array_nodes, 
				       rutils->array_or,
				       rutils->working_ca, 
				       MAX_LEN_MUT_ALLELE,
				       tmp1, tmp2, tmp3,
				       ignore_first, ignore_last, 
				       expected_covg, &m);


    }
  
  strbuf_free(tmp1);
  strbuf_free(tmp2);
  strbuf_free(tmp3);

}


//one gene per fasta, so if you have multiple reads,
//these are different exemplars, for divergent versions of the same gene
void load_antibiotic_gene_presence_info_on_sample(FILE* fp,
						  dBGraph* db_graph,
						  int (*file_reader)(FILE * fp, 
								     Sequence * seq, 
								     int max_read_length, 
								     boolean new_entry, 
								     boolean * full_entry),
						  AntibioticInfo* abi,
						  ReadingUtils* rutils,
						  GeneInfo* tmp_gi)



{
  reset_reading_utils(rutils);

  int num=1;

  while (num>0)
    {
      num = get_next_gene_info(fp,
			       db_graph,
			       tmp_gi,
			       rutils->seq,
			       rutils->kmer_window,
			       file_reader,
			       rutils->array_nodes,
			       rutils->array_or,
			       rutils->working_ca,
			       MAX_LEN_GENE);
      /*      printf("Percent >0 %d\n Median on nonzero %d\nMin %d\n, median %d\n",
	     tmp_gi->percent_nonzero,
	     tmp_gi->median_covg_on_nonzero_nodes,
	     tmp_gi->median_covg,
	     tmp_gi->min_covg);
      */

      if (tmp_gi->percent_nonzero>abi->genes[tmp_gi->name]->percent_nonzero)
	{
	  copy_gene_info(tmp_gi, abi->genes[tmp_gi->name]);
	}

    }
  

}



void load_antibiotic_mut_and_gene_info(dBGraph* db_graph,
				       int (*file_reader)(FILE * fp, 
							  Sequence * seq, 
							  int max_read_length, 
							  boolean new_entry, 
							  boolean * full_entry),
				       AntibioticInfo* abi,
				       ReadingUtils* rutils,
				       VarOnBackground* tmp_vob,
				       GeneInfo* tmp_gi,
				       int ignore_first, int ignore_last, int expected_covg,
				       StrBuf* install_dir)

{

  FILE* fp;

  if (abi->num_mutations>0)
    {
      fp = fopen(abi->m_fasta->buff, "r");
      if (fp==NULL)
	{
	  die("Cannot open %s - should be there as part of the install - did you run out of disk mid-install?\n",
	      abi->m_fasta->buff);
	}
      
      load_antibiotic_mutation_info_on_sample(fp,
					      db_graph,
					      file_reader,
					      abi,
					      rutils, 
					      tmp_vob,
					      ignore_first, ignore_last, expected_covg);
      fclose(fp);
    }
  if (abi->num_genes>0)
    {
      StrBuf* tmp = strbuf_new();
      int j;
      for (j=0; j<abi->num_genes; j++)
	{
	  strbuf_reset(tmp);
	  GenePresenceGene g = (GenePresenceGene) abi->which_genes[j];
	  map_gene_to_fasta(g, tmp, install_dir);
	  FILE* fp = fopen(tmp->buff, "r");
	  if (fp==NULL)
	    {
	      die("Unable to open %s, which should come as part of the install.\n", tmp->buff);
	    }
	  load_antibiotic_gene_presence_info_on_sample(fp,
						       db_graph,
						       file_reader,
						       abi,
						       rutils,
						       tmp_gi);
	  fclose(fp);
	}
      strbuf_free(tmp);

    }
}




InfectionType is_amikacin_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = amikacin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);
      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/amikacin.fa");
      abi->num_mutations = 4;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rrs_A1401X;
      int last_mut = rrs_G1484X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
	resistotype(abi->vars[i],
		    err_rate, 
		    db_graph->kmer_size, 
		    lambda_g, 
		    lambda_e, 
		    epsilon,
		    &best_model, 
		    MaxAPosteriori,
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
	  // Save the variant being called 
	  // add_called_variant_info_to_array(i, called_variants, best_model, I)
	  update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
	  return I;
	}
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_capreomycin_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = capreomycin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/capreomycin.fa");

      abi->num_mutations = 4;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rrs_A1401X;
      int last_mut = rrs_G1484X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
        err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_ethambutol_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = ethambutol;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/ethambutol.fa");

      abi->num_mutations = 1;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = embB_M306X;
      int last_mut = embB_M306X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
              err_rate, 
              db_graph->kmer_size, 
              lambda_g, 
              lambda_e, 
              epsilon,
              &best_model,
              MaxAPosteriori,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_isoniazid_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = isoniazid;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/isoniazid.fa");

      abi->num_mutations = 4;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = katG_S315X;
      int last_mut = fabG1_Gu17X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
        err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_kanamycin_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = kanamycin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/kanamycin.fa");

      abi->num_mutations = 4;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rrs_A1401X;
      int last_mut = rrs_G1484X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
        err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_pyrazinamide_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = pyrazinamide;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/pyrazinamide.fa");

      abi->num_mutations = 70;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = pncA_D49N;
      int last_mut = pncA_V21G;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
        err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_quinolones_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = quinolones;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/quinolones.fa");

      abi->num_mutations = 40;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = gyrA_H85X;
      int last_mut = gyrA_D94X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
        err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }

}
    InfectionType is_rifampicin_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = rifampicin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/rifampicin.fa");

      abi->num_mutations = 130;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rpoB_F425X;
      int last_mut = rpoB_L452X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
  double max_sus_conf=0;
  double min_conf=9999999999;
  boolean any_allele_non_null=false;
  boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
	{
	  continue;
	}
      any_allele_non_null=true;
      InfectionType I=
	resistotype(abi->vars[i],
		    err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
	  return I;
	}
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }

  }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }


}
    InfectionType is_streptomycin_susceptible(dBGraph* db_graph,
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
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate,
                 CalledVariant* called_variants,CalledGene* called_genes,
                 CmdLine* cmd_line )
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = streptomycin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/streptomycin.fa");

      abi->num_mutations = 99;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_vob,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rpsL_K43R;
      int last_mut = rrs_C517X;

      int i;
       //if you have any of these resistance alleles - call resistant
  Model best_model;
double max_sus_conf=0;
double min_conf=9999999999;
boolean any_allele_non_null=false;
boolean any_unsure_mixed_call=false;
  for (i=first_mut; i<=last_mut; i++)
    {

      if (both_alleles_null(abi->vars[i])==true)
  {
    continue;
  }
      any_allele_non_null=true;
      InfectionType I=
  resistotype(abi->vars[i],
        err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
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
    return I;
  }
    if (I == Unsure && best_model.type == MixedInfection)
  {
    any_unsure_mixed_call = true;
     update_called_variants(called_variants,i,abi->vars[i], best_model.conf);
  }
    }
  if (
      (any_allele_non_null==false) || (any_unsure_mixed_call)
      )

    {
      return Susceptible;
    }
  else if (max_sus_conf>MIN_CONFIDENCE_S)
    {
      return Susceptible;
    }
  else
    {
      return Unsure;
    }


}

void print_antibiotic_susceptibility(dBGraph* db_graph,
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
                      int ignore_first, 
                      int ignore_last, 
                      int expected_covg,
                      double lambda_g, 
                      double lambda_e, 
                      double err_rate,
                      CalledVariant* called_variants,
                      CalledGene* called_genes,
                      CmdLine* cmd_line),
          StrBuf* tmpbuf,
          StrBuf* install_dir,
          int ignore_first, int ignore_last,
          int expected_covg,
          double lambda_g, double lambda_e, double err_rate,
          CmdLine* cmd_line,
          boolean output_last,//for JSON
          CalledVariant* called_variants,
          CalledGene* called_genes
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
        ignore_first, 
        ignore_last, 
        expected_covg,
        lambda_g,
        lambda_e,
        err_rate,
        called_variants,
        called_genes,
        cmd_line);

  
  map_antibiotic_enum_to_str(abi->ab, tmpbuf);
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
      else if (  suc==Resistant )
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
