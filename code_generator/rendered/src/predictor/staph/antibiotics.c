/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *  antibiotics.c 
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

   switch (gene) 
   {
    
    case vgbA : return "Clindamycin";
    
    case IsaB : return "Clindamycin";
    
    case blaZ : return "Penicillin";
    
    case lnuA : return "Clindamycin";
    
    case lnuB : return "Clindamycin";
    
    case ermB : return "Erythromycin";
    
    case ermC : return "Erythromycin";
    
    case ermA : return "Erythromycin";
    
    case ermY : return "Erythromycin";
    
    case mupB : return "Mupirocin";
    
    case mupA : return "Mupirocin";
    
    case ermT : return "Erythromycin";
    
    case dfrA : return "Trimethoprim";
    
    case vgaALC : return "Clindamycin";
    
    case dfrC : return "Trimethoprim";
    
    case dfrD : return "Trimethoprim";
    
    case dfrG : return "Trimethoprim";
    
    case mecC : return "Methicillin";
    
    case mecA : return "Methicillin";
    
    case dfrK : return "Trimethoprim";
    
    case vgaB : return "Clindamycin";
    
    case vgaA : return "Clindamycin";
    
    case tetK : return "Tetracycline";
    
    case tetM : return "Tetracycline";
    
    case tetL : return "Tetracycline";
    
    case tetO : return "Tetracycline";
    
    case msrA : return "Erythromycin";
    
    case vanA : return "Vancomycin";
    
    case vanC : return "Vancomycin";
    
    case vanB : return "Vancomycin";
    
    case fusB : return "FusidicAcid";
    
    case fusC : return "FusidicAcid";
    
    case str : return "Gentamicin";
    
    case aacAaphD : return "Gentamicin";
    
    case PVL : return "";
    
    case unspecified_gpg : return "unknown";    
   }
	   
   return "unknown";
}
void map_antibiotic_enum_to_str(Antibiotic ab, StrBuf* name)
{
  if (ab==NoDrug)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "NoDrug");
    }
  
  else if (ab==Erythromycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Erythromycin");
    }
  
  else if (ab==Clindamycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Clindamycin");
    }
  
  else if (ab==Penicillin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Penicillin");
    }
  
  else if (ab==Mupirocin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Mupirocin");
    }
  
  else if (ab==Trimethoprim)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Trimethoprim");
    }
  
  else if (ab==Methicillin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Methicillin");
    }
  
  else if (ab==Tetracycline)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Tetracycline");
    }
  
  else if (ab==Vancomycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Vancomycin");
    }
  
  else if (ab==FusidicAcid)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "FusidicAcid");
    }
  
  else if (ab==Gentamicin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Gentamicin");
    }
  
  else if (ab==Rifampicin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Rifampicin");
    }
  
  else if (ab==Ciprofloxacin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Ciprofloxacin");
    }
  
  else
    {
      die("Impossible - compiler should not allow this\n");
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
  // TODO fix issue causing bug here. If we uncomment this with no mutations we get a seg fault
  // for (i=0; i<NUM_KNOWN_MUTATIONS; i++)
  //   {
  //     reset_var_on_background(abi->vars[i]->vob_best_sus);
  //     reset_var_on_background(abi->vars[i]->vob_best_res);
  //   }
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
					      int ignore_first, int ignore_last)
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
				       &m);


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
				       int ignore_first, 
				       int ignore_last, 
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
					      ignore_first, ignore_last);
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




void update_infection_type(InfectionType* I_new, InfectionType* I_permenant){
  if ( (*I_permenant==Unsure) || (*I_permenant==Susceptible) ){
    *I_permenant = *I_new;
  }
}


InfectionType is_erythromycin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
   *any_erm_present=false;


  //setup antibiotic info object
  abi->ab = Erythromycin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/erythromycin.fa");
  
    abi->which_genes[0]=ermB;
  
    abi->which_genes[1]=ermC;
  
    abi->which_genes[2]=ermA;
  
    abi->which_genes[3]=ermY;
  
    abi->which_genes[4]=ermT;
  
    abi->which_genes[5]=msrA;
  
  abi->num_genes=6;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<6; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN_ERY 
    ,&genotyped_present

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
         
          if ( (abi->which_genes[i] == ermA) || (abi->which_genes[i] == ermB) || (abi->which_genes[i] == ermC) || (abi->which_genes[i] == ermY) || (abi->which_genes[i] == ermT) )
          {
            *any_erm_present=true;
          }
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

}


InfectionType is_clindamycin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Clindamycin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/clindamycin.fa");
  
    abi->which_genes[0]=vgbA;
  
    abi->which_genes[1]=IsaB;
  
    abi->which_genes[2]=lnuA;
  
    abi->which_genes[3]=lnuB;
  
    abi->which_genes[4]=vgaALC;
  
    abi->which_genes[5]=vgaB;
  
    abi->which_genes[6]=vgaA;
  
  abi->num_genes=7;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<7; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

}


InfectionType is_penicillin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Penicillin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/penicillin.fa");
  
    abi->which_genes[0]=blaZ;
  
  abi->num_genes=1;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;


genotyped_present = false;
InfectionType I= resistotype_gene(abi->genes[blaZ], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
         &best_model, MaxAPosteriori,
          MIN_PERC_COVG_BLAZ ,

          MIN_GENE_CN_PEN
		,&genotyped_present
         );
  if ( genotyped_present ||  cmd_line->verbose) {
    update_called_genes(called_genes, blaZ, abi->genes[blaZ], best_model.conf );
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
  


  
  
   
  return I_permanent;
  

}


InfectionType is_mupirocin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Mupirocin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/mupirocin.fa");
  
    abi->which_genes[0]=mupB;
  
    abi->which_genes[1]=mupA;
  
  abi->num_genes=2;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<2; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN_MUP 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

}


InfectionType is_trimethoprim_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Trimethoprim;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/trimethoprim.fa");
  
    abi->which_genes[0]=dfrA;
  
    abi->which_genes[1]=dfrC;
  
    abi->which_genes[2]=dfrD;
  
    abi->which_genes[3]=dfrG;
  
    abi->which_genes[4]=dfrK;
  
  abi->num_genes=5;
  abi->num_mutations = 8;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;

  int first_mut = dfrB_F99I;
  int last_mut = dfrB_N60I;
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
      genotyped_present = false;
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
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
    }



for (i=0; i<5; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_DFRK ,
    
     MIN_GENE_CN_PEN 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
  
  if( (I_permanent==Resistant) || (I_permanent==MixedInfection) ) {
    return I_permanent;
  }
  else{
    if (any_allele_non_null==false)
      {
        return Unsure;
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
  

}


InfectionType is_methicillin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Methicillin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/methicillin.fa");
  
    abi->which_genes[0]=mecC;
  
    abi->which_genes[1]=mecA;
  
  abi->num_genes=2;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<2; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN_MEC 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

}


InfectionType is_tetracycline_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Tetracycline;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/tetracycline.fa");
  
    abi->which_genes[0]=tetK;
  
    abi->which_genes[1]=tetM;
  
    abi->which_genes[2]=tetL;
  
    abi->which_genes[3]=tetO;
  
  abi->num_genes=4;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<4; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN_TET 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

}


InfectionType is_vancomycin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Vancomycin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/vancomycin.fa");
  
    abi->which_genes[0]=vanA;
  
    abi->which_genes[1]=vanC;
  
    abi->which_genes[2]=vanB;
  
  abi->num_genes=3;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<3; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

}


InfectionType is_fusidicacid_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = FusidicAcid;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/fusidicacid.fa");
  
    abi->which_genes[0]=fusB;
  
    abi->which_genes[1]=fusC;
  
  abi->num_genes=2;
  abi->num_mutations = 43;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;

  int first_mut = fusA_A655P;
  int last_mut = fusA_E468V;
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
      genotyped_present = false;
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
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
    }



for (i=0; i<2; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN_FUS 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  



  genotyped_present = false;
  InfectionType I_f652s=
  resistotype(abi->vars[fusA_F652S],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);

  genotyped_present = false;
  InfectionType I_y654n=
    resistotype(abi->vars[fusA_Y654N],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);
if (I_f652s==Resistant && I_y654n==Resistant)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_F652S], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_Y654N], best_model.conf);
    update_infection_type(Resistant,&I_permanent);  
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_F652S], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_Y654N], best_model.conf);    
  }




  genotyped_present = false;
  InfectionType I_t326i=
  resistotype(abi->vars[fusA_T326I],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg,contaminiation_covg, 
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);

  genotyped_present = false;
  InfectionType I_e468v=
    resistotype(abi->vars[fusA_E468V],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);


if (I_t326i==Resistant && I_e468v==Resistant)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_T326I],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E468V],best_model.conf);
    update_infection_type(Resistant,&I_permanent);  
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_T326I],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E468V],best_model.conf);
  }

  



  genotyped_present = false;
  InfectionType I_l461f=
  resistotype(abi->vars[fusA_L461F],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);

  genotyped_present = false;
  InfectionType I_a376v=
    resistotype(abi->vars[fusA_A376V],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);

  genotyped_present = false;
  InfectionType I_a655p=
    resistotype(abi->vars[fusA_A655P],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);

  genotyped_present = false;
  InfectionType I_d463g=
    resistotype(abi->vars[fusA_D463G],
		err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
		&best_model, MaxAPosteriori,
		cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);
  
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
    update_infection_type(Resistant,&I_permanent);  
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_A376V],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_A655P],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_D463G],best_model.conf);
  }

genotyped_present = false;
InfectionType I_e444v=Susceptible;
resistotype(abi->vars[fusA_E444V],
	    err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
	    &best_model, MaxAPosteriori,
	    cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);
if ((I_l461f==Resistant)
       &&
    (I_e444v==Resistant) )
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E444V],best_model.conf); 
    update_infection_type(Resistant,&I_permanent);     
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[fusA_L461F],best_model.conf);
    update_called_variants(called_variants,i,abi->vars[fusA_E444V],best_model.conf); 
  }  



  
  
  
  if( (I_permanent==Resistant) || (I_permanent==MixedInfection) ) {
    return I_permanent;
  }
  else{
    if (any_allele_non_null==false)
      {
        return Unsure;
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
  

}


InfectionType is_gentamicin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Gentamicin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/gentamicin.fa");
  
    abi->which_genes[0]=str;
  
    abi->which_genes[1]=aacAaphD;
  
  abi->num_genes=2;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;



for (i=0; i<2; i++)
    {
      genotyped_present = false;
      InfectionType I =
	     resistotype_gene(abi->genes[abi->which_genes[i]], 
			 err_rate, db_graph->kmer_size, 
			 lambda_g, lambda_e, epsilon, expected_covg, contaminiation_covg,
			 &best_model, MaxAPosteriori,
			  MIN_PERC_COVG_STANDARD ,
    
     MIN_GENE_CN_GEN 
    ,&genotyped_present

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
        
      }
        if ( genotyped_present ||  cmd_line->verbose) {
          update_called_genes(called_genes, abi->which_genes[i], abi->genes[abi->which_genes[i]], best_model.conf );
        }      
        update_infection_type(&I,&I_permanent);
    }
  


  
  
   
  return I_permanent;
  

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
            int ignore_first, int ignore_last, SpeciesInfo* species_info,
            double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Rifampicin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/rifampicin.fa");
  
  abi->num_genes=0;
  abi->num_mutations = 20;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;

  int first_mut = rpoB_A477D;
  int last_mut = rpoB_N474K;
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
      genotyped_present = false;
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
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
    }

  


  genotyped_present = false;
  InfectionType I_m470t=
  resistotype(abi->vars[rpoB_M470T],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops,
    &genotyped_present);

if ( (I_m470t==Susceptible) && (best_model.conf>max_sus_conf) )
  {
    max_sus_conf = best_model.conf;
  }
if (best_model.conf<min_conf)
  {
    min_conf = best_model.conf;
  }

  genotyped_present = false;
  InfectionType I_d471g=
  resistotype(abi->vars[rpoB_D471G],
	      err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
	      &best_model, MaxAPosteriori,
	      cmd_line->min_frac_to_detect_minor_pops,
        &genotyped_present);

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
    update_infection_type(Resistant,&I_permanent);   //ignoring mixed infections for epistatic case
  }
  else if (cmd_line->verbose)
  {
    update_called_variants(called_variants,i,abi->vars[rpoB_D471G], best_model.conf);
    update_called_variants(called_variants,i,abi->vars[rpoB_M470T], best_model.conf);
  }    
  



  
  
  
  if( (I_permanent==Resistant) || (I_permanent==MixedInfection) ) {
    return I_permanent;
  }
  else{
    if (any_allele_non_null==false)
      {
        return Unsure;
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
  

}


InfectionType is_ciprofloxacin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
            )

{
  
    int expected_covg = species_info->species_covg_info->median_coverage[Saureus];
  
  int contaminiation_covg = get_contamination_covg(species_info);

  InfectionType I_permanent = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Ciprofloxacin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/ciprofloxacin.fa");
  
  abi->num_genes=0;
  abi->num_mutations = 6;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
            file_reader,
            abi,
            rutils,
            tmp_vob,
            tmp_gi,
            ignore_first, ignore_last,
            install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  boolean genotyped_present = false;

  int first_mut = gyrA_E88K;
  int last_mut = grlA_S80Y;
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
      genotyped_present = false;
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,expected_covg, contaminiation_covg,
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
    }

  


  
  
  
  if( (I_permanent==Resistant) || (I_permanent==MixedInfection) ) {
    return I_permanent;
  }
  else{
    if (any_allele_non_null==false)
      {
        return Unsure;
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
                    SpeciesInfo* species_info,
      							double lambda_g,
                    double lambda_e,
                    double err_rate,
                    CalledVariant* called_variants,
                    CalledGene* called_genes,
                    CmdLine* cmd_line),
					StrBuf* tmpbuf,
					StrBuf* install_dir,
					int ignore_first, int ignore_last,
					SpeciesInfo* species_info,
					double lambda_g, double lambda_e, double err_rate,
          CmdLine* cmd_line,
					boolean output_last,//for JSON,
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
	      species_info,
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
      else if ( suc==Resistant )
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





///virulence






    
	Troolean is_pvl_positive(dBGraph* db_graph,
			   int (*file_reader)(FILE * fp, 
					      Sequence * seq, 
					      int max_read_length, 
					      boolean new_entry, 
					      boolean * full_entry),
			ReadingUtils* rutils,
			GeneInfo* tmp_gi,
			StrBuf* install_dir)
			   

{
  StrBuf* fa = strbuf_create(install_dir->buff);
  strbuf_append_str(fa, "data/staph/virulence/PVL.fa");

  FILE* fp = fopen(fa->buff, "r");
  if (fp==NULL)
    {
      die("Cannot open %s\n", fa->buff);
    }
  int num=1;
  boolean is_pos=false;
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
      if (tmp_gi->percent_nonzero > MIN_PERC_COVG_VIRULENCE)
	{
	  is_pos=true;
	}
    }
  fclose(fp);
  strbuf_free(fa);
  return is_pos;

}
	void print_pvl_presence(dBGraph* db_graph,
			int (*file_reader)(FILE * fp, 
					   Sequence * seq, 
					   int max_read_length, 
					   boolean new_entry, 
					   boolean * full_entry),
			ReadingUtils* rutils,
			GeneInfo* tmp_gi,
			Troolean (*func)(dBGraph* db_graph,
					int (*file_reader)(FILE * fp, 
							   Sequence * seq, 
							   int max_read_length, 
							   boolean new_entry, 
							   boolean * full_entry),
					ReadingUtils* rutils,
					GeneInfo* tmp_gi,
					StrBuf* install_dir),
			StrBuf* install_dir, OutputFormat format)
{

  Troolean result = is_pvl_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("PVL\t");
      if (result==true)
	{
	  printf("positive\n");
	}
      else
	{
	  printf("negative\n");
	}
    }
  else
    {
      
      if (result==true)
	{
	  print_json_item("PVL", "positive", true  );
	}
      else
	{
	  print_json_item("PVL", "negative", true );
	}
      
    }
}




