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
    
    case aacAaphD : return "Gentamicin";
    
    case IsaB : return "Clindamycin";
    
    case aadEant6Ia : return "Gentamicin";
    
    case aadDaph4Ia : return "Gentamicin";
    
    case qacB : return "Biocides";
    
    case qacA : return "Biocides";
    
    case cfr : return "Chloramphenicol,Linezolid";
    
    case blaZ : return "Penicillin";
    
    case lnuA : return "Clindamycin";
    
    case lnuB : return "Clindamycin";
    
    case ermB : return "Erythromycin";
    
    case ermC : return "Erythromycin";
    
    case ermA : return "Erythromycin";
    
    case ermY : return "Erythromycin";
    
    case mupB : return "Mupirocin";
    
    case mupA : return "Mupirocin";
    
    case aph2Ic : return "Gentamicin";
    
    case ant9Ia : return "Gentamicin";
    
    case ermT : return "Erythromycin";
    
    case ant9Ib : return "Spectinomycin";
    
    case aphA3aph3III : return "Gentamicin";
    
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
    
    case mphC : return "Erythromycin";
    
    case sat4 : return "Streptothricin";
    
    case cat : return "Chloramphenicol";
    
    case vanA : return "Vancomycin";
    
    case vanC : return "Vancomycin";
    
    case vanB : return "Vancomycin";
    
    case fusB : return "FusidicAcid";
    
    case fusC : return "FusidicAcid";
    
    case str : return "Gentamicin";
    
    case qacCsmr : return "Biocides";
    
    case arcA : return "";
    
    case arcB : return "";
    
    case arcC : return "";
    
    case arcD : return "";
    
    case ccrA : return "";
    
    case ccrB : return "";
    
    case ccrCa : return "";
    
    case ccrCb : return "";
    
    case ccrCc : return "";
    
    case chp : return "";
    
    case eta : return "";
    
    case etb : return "";
    
    case etd : return "";
    
    case luk : return "";
    
    case lukM : return "";
    
    case lukMF : return "";
    
    case lukPVF : return "";
    
    case lukPVS : return "";
    
    case sak : return "";
    
    case sasX : return "";
    
    case scn : return "";
    
    case sea : return "";
    
    case seb : return "";
    
    case sec : return "";
    
    case sed : return "";
    
    case see : return "";
    
    case seg : return "";
    
    case seh : return "";
    
    case sei : return "";
    
    case sej : return "";
    
    case selR : return "";
    
    case sep : return "";
    
    case seu : return "";
    
    case tsst1 : return "";
    
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
  
  else if (ab==Gentamicin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Gentamicin");
    }
  
  else if (ab==Biocides)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Biocides");
    }
  
  else if (ab==Chloramphenicol)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Chloramphenicol");
    }
  
  else if (ab==Linezolid)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Linezolid");
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
  
  else if (ab==Spectinomycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Spectinomycin");
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
  
  else if (ab==Streptothricin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Streptothricin");
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
				       int ignore_first, 
				       int ignore_last, 
				       int expected_covg,
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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
             boolean* any_erm_present,
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
  
    abi->which_genes[6]=mphC;
  
  abi->num_genes=7;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<7; i++)
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
         
          if ( (abi->which_genes[i] == ermA) || (abi->which_genes[i] == ermB) || (abi->which_genes[i] == ermC) || (abi->which_genes[i] == ermY) || (abi->which_genes[i] == ermT) )
          {
            *any_erm_present=true;
          }
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<7; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Gentamicin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/gentamicin.fa");
  
    abi->which_genes[0]=aacAaphD;
  
    abi->which_genes[1]=aadEant6Ia;
  
    abi->which_genes[2]=aadDaph4Ia;
  
    abi->which_genes[3]=aph2Ic;
  
    abi->which_genes[4]=ant9Ia;
  
    abi->which_genes[5]=aphA3aph3III;
  
    abi->which_genes[6]=str;
  
  abi->num_genes=7;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<7; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

}


InfectionType is_biocides_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Biocides;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/biocides.fa");
  
    abi->which_genes[0]=qacB;
  
    abi->which_genes[1]=qacA;
  
    abi->which_genes[2]=qacCsmr;
  
  abi->num_genes=3;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<3; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

}


InfectionType is_chloramphenicol_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Chloramphenicol;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/chloramphenicol.fa");
  
    abi->which_genes[0]=cfr;
  
    abi->which_genes[1]=cat;
  
  abi->num_genes=2;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<2; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

}


InfectionType is_linezolid_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Linezolid;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/linezolid.fa");
  
    abi->which_genes[0]=cfr;
  
  abi->num_genes=1;
  abi->num_mutations = 1;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;

  int mut = twentythreeS_G2576T;
  i = mut;
  boolean any_allele_non_null=false;
  if (both_alleles_null(abi->vars[i])==true)
	{
	  any_allele_non_null=true;
	}
  I=resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
    lambda_g, lambda_e, epsilon,
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
	}
  update_infection_type(&I,&I_permenant);


I= resistotype_gene(abi->genes[cfr], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon,expected_covg,
         &best_model, MaxAPosteriori,
         MIN_PERC_COVG_BLAZ);
  if ( (I==Resistant) || (I==MixedInfection) ) {
    update_called_genes(called_genes, cfr, abi->genes[cfr], best_model.conf );
  }
  update_infection_type(&I,&I_permenant);
  


  
  
  
  if( (I_permenant==Resistant) || (I_permenant==MixedInfection) ) {
    return I_permenant;
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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;


I= resistotype_gene(abi->genes[blaZ], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon,expected_covg,
         &best_model, MaxAPosteriori,
         MIN_PERC_COVG_BLAZ);
  if ( (I==Resistant) || (I==MixedInfection) ) {
    update_called_genes(called_genes, blaZ, abi->genes[blaZ], best_model.conf );
  }
  update_infection_type(&I,&I_permenant);
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<2; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

}


InfectionType is_spectinomycin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Spectinomycin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/spectinomycin.fa");
  
    abi->which_genes[0]=ant9Ib;
  
  abi->num_genes=1;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;


I= resistotype_gene(abi->genes[ant9Ib], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon,expected_covg,
         &best_model, MaxAPosteriori,
         MIN_PERC_COVG_BLAZ);
  if ( (I==Resistant) || (I==MixedInfection) ) {
    update_called_genes(called_genes, ant9Ib, abi->genes[ant9Ib], best_model.conf );
  }
  update_infection_type(&I,&I_permenant);
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;

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
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,
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
    	}
      update_infection_type(&I,&I_permenant);
    }



for (i=0; i<5; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
  
  if( (I_permenant==Resistant) || (I_permenant==MixedInfection) ) {
    return I_permenant;
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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<2; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<4; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

}


InfectionType is_streptothricin_susceptible(dBGraph* db_graph,
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
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
  reset_antibiotic_info(abi);
  


  //setup antibiotic info object
  abi->ab = Streptothricin;
  strbuf_append_str(abi->m_fasta, install_dir->buff);
  strbuf_append_str(abi->m_fasta, "data/staph/antibiotics/streptothricin.fa");
  
    abi->which_genes[0]=sat4;
  
  abi->num_genes=1;
  abi->num_mutations = 0;

  double epsilon = pow(1-err_rate, db_graph->kmer_size);
  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_vob,
				    tmp_gi,
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;


I= resistotype_gene(abi->genes[sat4], err_rate, db_graph->kmer_size, 
         lambda_g, lambda_e, epsilon,expected_covg,
         &best_model, MaxAPosteriori,
         MIN_PERC_COVG_BLAZ);
  if ( (I==Resistant) || (I==MixedInfection) ) {
    update_called_genes(called_genes, sat4, abi->genes[sat4], best_model.conf );
  }
  update_infection_type(&I,&I_permenant);
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;



for (i=0; i<3; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


  
  
   
  return I_permenant;
  

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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;

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
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,
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
    	}
      update_infection_type(&I,&I_permenant);
    }



for (i=0; i<2; i++)
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
        
          update_called_genes(called_genes,  abi->which_genes[i] , abi->genes[abi->which_genes[i]],best_model.conf );
        }
        update_infection_type(&I,&I_permenant);
    }
  


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



  
  
  
  if( (I_permenant==Resistant) || (I_permenant==MixedInfection) ) {
    return I_permenant;
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
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;

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
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,
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
    	}
      update_infection_type(&I,&I_permenant);
    }

  


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



  
  
  
  if( (I_permenant==Resistant) || (I_permenant==MixedInfection) ) {
    return I_permenant;
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
				    int ignore_first, int ignore_last, int expected_covg,
				    double lambda_g, double lambda_e, double err_rate,
            
             CalledVariant* called_variants,CalledGene* called_genes,
             CmdLine* cmd_line
				    )

{
  InfectionType I_permenant = Unsure;
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
				    ignore_first, ignore_last, expected_covg,
				    install_dir);
  double max_sus_conf=0;
  double min_conf=9999999;  
  int i;
  Model best_model;
  InfectionType I;

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
      InfectionType I=
	    resistotype(abi->vars[i], err_rate, db_graph->kmer_size, 
		    lambda_g, lambda_e, epsilon,
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
    	}
      update_infection_type(&I,&I_permenant);
    }

  


  
  
  
  if( (I_permenant==Resistant) || (I_permenant==MixedInfection) ) {
    return I_permenant;
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
							  int ignore_first, int ignore_last, int expected_covg,
							  double lambda_g, double lambda_e, double err_rate, 
							  boolean* any_erm_present,
                CalledVariant* called_variants,CalledGene* called_genes,
                CmdLine* cmd_line),
					  StrBuf* tmpbuf,
					  StrBuf* install_dir,
					  int ignore_first, int ignore_last, int expected_covg,
					  double lambda_g, double lambda_e, double err_rate, CmdLine* cmd_line, boolean output_last,//for JSON 
					  boolean* any_erm_present,
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
	      ignore_first, ignore_last, expected_covg,
	      lambda_g,
	      lambda_e,
	      err_rate,
	      any_erm_present,
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
							 int ignore_first, int ignore_last, int expected_covg,
							 double lambda_g, double lambda_e, double err_rate,
               CalledVariant* called_variants,CalledGene* called_genes,
               CmdLine* cmd_line),
					 StrBuf* tmpbuf,
					 boolean any_erm_present,
					 StrBuf* install_dir,
					 int ignore_first, int ignore_last, int expected_covg,
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
	      ignore_first, ignore_last, expected_covg,
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
      //the ordering of these if's matters
      if (suc==MixedInfection)
	{
	  printf("r(constitutive)\n");
	}
      else if (suc==Resistant)
	{
	  printf("R(constitutive)\n");
	}
      else if ( (suc==Susceptible) && (any_erm_present==true) )
	{
	  printf("R(inducible)\n");
	}
      else if (suc==Susceptible)
	{
	  printf("S\n");
	}
      else
	{
	  printf("N\n");
	}

    }
  else
    {
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
	  print_json_item(tmpbuf->buff, "R(inducible)", output_last);
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

}




///virulence






    
	Troolean is_arca_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/arcA.fa");

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
	void print_arca_presence(dBGraph* db_graph,
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

  Troolean result = is_arca_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("arcA\t");
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
	  print_json_item("arcA", "positive",  false  );
	}
      else
	{
	  print_json_item("arcA", "negative",  false );
	}
      
    }
}



    
	Troolean is_arcb_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/arcB.fa");

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
	void print_arcb_presence(dBGraph* db_graph,
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

  Troolean result = is_arcb_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("arcB\t");
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
	  print_json_item("arcB", "positive",  false  );
	}
      else
	{
	  print_json_item("arcB", "negative",  false );
	}
      
    }
}



    
	Troolean is_arcc_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/arcC.fa");

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
	void print_arcc_presence(dBGraph* db_graph,
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

  Troolean result = is_arcc_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("arcC\t");
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
	  print_json_item("arcC", "positive",  false  );
	}
      else
	{
	  print_json_item("arcC", "negative",  false );
	}
      
    }
}



    
	Troolean is_arcd_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/arcD.fa");

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
	void print_arcd_presence(dBGraph* db_graph,
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

  Troolean result = is_arcd_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("arcD\t");
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
	  print_json_item("arcD", "positive",  false  );
	}
      else
	{
	  print_json_item("arcD", "negative",  false );
	}
      
    }
}



    
	Troolean is_ccra_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/ccrA.fa");

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
	void print_ccra_presence(dBGraph* db_graph,
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

  Troolean result = is_ccra_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("ccrA\t");
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
	  print_json_item("ccrA", "positive",  false  );
	}
      else
	{
	  print_json_item("ccrA", "negative",  false );
	}
      
    }
}



    
	Troolean is_ccrb_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/ccrB.fa");

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
	void print_ccrb_presence(dBGraph* db_graph,
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

  Troolean result = is_ccrb_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("ccrB\t");
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
	  print_json_item("ccrB", "positive",  false  );
	}
      else
	{
	  print_json_item("ccrB", "negative",  false );
	}
      
    }
}



    
	Troolean is_ccrca_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/ccrCa.fa");

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
	void print_ccrca_presence(dBGraph* db_graph,
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

  Troolean result = is_ccrca_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("ccrCa\t");
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
	  print_json_item("ccrCa", "positive",  false  );
	}
      else
	{
	  print_json_item("ccrCa", "negative",  false );
	}
      
    }
}



    
	Troolean is_ccrcb_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/ccrCb.fa");

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
	void print_ccrcb_presence(dBGraph* db_graph,
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

  Troolean result = is_ccrcb_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("ccrCb\t");
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
	  print_json_item("ccrCb", "positive",  false  );
	}
      else
	{
	  print_json_item("ccrCb", "negative",  false );
	}
      
    }
}



    
	Troolean is_ccrcc_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/ccrCc.fa");

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
	void print_ccrcc_presence(dBGraph* db_graph,
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

  Troolean result = is_ccrcc_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("ccrCc\t");
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
	  print_json_item("ccrCc", "positive",  false  );
	}
      else
	{
	  print_json_item("ccrCc", "negative",  false );
	}
      
    }
}



    
	Troolean is_chp_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/chp.fa");

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
	void print_chp_presence(dBGraph* db_graph,
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

  Troolean result = is_chp_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("chp\t");
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
	  print_json_item("chp", "positive",  false  );
	}
      else
	{
	  print_json_item("chp", "negative",  false );
	}
      
    }
}



    
	Troolean is_eta_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/eta.fa");

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
	void print_eta_presence(dBGraph* db_graph,
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

  Troolean result = is_eta_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("eta\t");
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
	  print_json_item("eta", "positive",  false  );
	}
      else
	{
	  print_json_item("eta", "negative",  false );
	}
      
    }
}



    
	Troolean is_etb_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/etb.fa");

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
	void print_etb_presence(dBGraph* db_graph,
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

  Troolean result = is_etb_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("etb\t");
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
	  print_json_item("etb", "positive",  false  );
	}
      else
	{
	  print_json_item("etb", "negative",  false );
	}
      
    }
}



    
	Troolean is_etd_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/etd.fa");

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
	void print_etd_presence(dBGraph* db_graph,
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

  Troolean result = is_etd_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("etd\t");
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
	  print_json_item("etd", "positive",  false  );
	}
      else
	{
	  print_json_item("etd", "negative",  false );
	}
      
    }
}



    
	Troolean is_luk_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/luk.fa");

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
	void print_luk_presence(dBGraph* db_graph,
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

  Troolean result = is_luk_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("luk\t");
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
	  print_json_item("luk", "positive",  false  );
	}
      else
	{
	  print_json_item("luk", "negative",  false );
	}
      
    }
}



    
	Troolean is_lukm_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/lukM.fa");

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
	void print_lukm_presence(dBGraph* db_graph,
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

  Troolean result = is_lukm_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("lukM\t");
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
	  print_json_item("lukM", "positive",  false  );
	}
      else
	{
	  print_json_item("lukM", "negative",  false );
	}
      
    }
}



    
	Troolean is_lukmf_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/lukMF.fa");

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
	void print_lukmf_presence(dBGraph* db_graph,
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

  Troolean result = is_lukmf_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("lukMF\t");
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
	  print_json_item("lukMF", "positive",  false  );
	}
      else
	{
	  print_json_item("lukMF", "negative",  false );
	}
      
    }
}



    
	Troolean is_lukpvf_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/lukPVF.fa");

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
	void print_lukpvf_presence(dBGraph* db_graph,
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

  Troolean result = is_lukpvf_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("lukPVF\t");
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
	  print_json_item("lukPVF", "positive",  false  );
	}
      else
	{
	  print_json_item("lukPVF", "negative",  false );
	}
      
    }
}



    
	Troolean is_lukpvs_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/lukPVS.fa");

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
	void print_lukpvs_presence(dBGraph* db_graph,
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

  Troolean result = is_lukpvs_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("lukPVS\t");
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
	  print_json_item("lukPVS", "positive",  false  );
	}
      else
	{
	  print_json_item("lukPVS", "negative",  false );
	}
      
    }
}



    
	Troolean is_sak_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sak.fa");

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
	void print_sak_presence(dBGraph* db_graph,
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

  Troolean result = is_sak_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sak\t");
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
	  print_json_item("sak", "positive",  false  );
	}
      else
	{
	  print_json_item("sak", "negative",  false );
	}
      
    }
}



    
	Troolean is_sasx_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sasX.fa");

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
	void print_sasx_presence(dBGraph* db_graph,
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

  Troolean result = is_sasx_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sasX\t");
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
	  print_json_item("sasX", "positive",  false  );
	}
      else
	{
	  print_json_item("sasX", "negative",  false );
	}
      
    }
}



    
	Troolean is_scn_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/scn.fa");

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
	void print_scn_presence(dBGraph* db_graph,
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

  Troolean result = is_scn_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("scn\t");
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
	  print_json_item("scn", "positive",  false  );
	}
      else
	{
	  print_json_item("scn", "negative",  false );
	}
      
    }
}



    
	Troolean is_sea_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sea.fa");

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
	void print_sea_presence(dBGraph* db_graph,
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

  Troolean result = is_sea_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sea\t");
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
	  print_json_item("sea", "positive",  false  );
	}
      else
	{
	  print_json_item("sea", "negative",  false );
	}
      
    }
}



    
	Troolean is_seb_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/seb.fa");

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
	void print_seb_presence(dBGraph* db_graph,
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

  Troolean result = is_seb_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("seb\t");
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
	  print_json_item("seb", "positive",  false  );
	}
      else
	{
	  print_json_item("seb", "negative",  false );
	}
      
    }
}



    
	Troolean is_sec_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sec.fa");

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
	void print_sec_presence(dBGraph* db_graph,
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

  Troolean result = is_sec_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sec\t");
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
	  print_json_item("sec", "positive",  false  );
	}
      else
	{
	  print_json_item("sec", "negative",  false );
	}
      
    }
}



    
	Troolean is_sed_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sed.fa");

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
	void print_sed_presence(dBGraph* db_graph,
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

  Troolean result = is_sed_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sed\t");
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
	  print_json_item("sed", "positive",  false  );
	}
      else
	{
	  print_json_item("sed", "negative",  false );
	}
      
    }
}



    
	Troolean is_see_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/see.fa");

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
	void print_see_presence(dBGraph* db_graph,
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

  Troolean result = is_see_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("see\t");
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
	  print_json_item("see", "positive",  false  );
	}
      else
	{
	  print_json_item("see", "negative",  false );
	}
      
    }
}



    
	Troolean is_seg_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/seg.fa");

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
	void print_seg_presence(dBGraph* db_graph,
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

  Troolean result = is_seg_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("seg\t");
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
	  print_json_item("seg", "positive",  false  );
	}
      else
	{
	  print_json_item("seg", "negative",  false );
	}
      
    }
}



    
	Troolean is_seh_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/seh.fa");

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
	void print_seh_presence(dBGraph* db_graph,
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

  Troolean result = is_seh_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("seh\t");
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
	  print_json_item("seh", "positive",  false  );
	}
      else
	{
	  print_json_item("seh", "negative",  false );
	}
      
    }
}



    
	Troolean is_sei_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sei.fa");

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
	void print_sei_presence(dBGraph* db_graph,
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

  Troolean result = is_sei_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sei\t");
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
	  print_json_item("sei", "positive",  false  );
	}
      else
	{
	  print_json_item("sei", "negative",  false );
	}
      
    }
}



    
	Troolean is_sej_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sej.fa");

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
	void print_sej_presence(dBGraph* db_graph,
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

  Troolean result = is_sej_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sej\t");
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
	  print_json_item("sej", "positive",  false  );
	}
      else
	{
	  print_json_item("sej", "negative",  false );
	}
      
    }
}



    
	Troolean is_selr_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/selR.fa");

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
	void print_selr_presence(dBGraph* db_graph,
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

  Troolean result = is_selr_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("selR\t");
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
	  print_json_item("selR", "positive",  false  );
	}
      else
	{
	  print_json_item("selR", "negative",  false );
	}
      
    }
}



    
	Troolean is_sep_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/sep.fa");

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
	void print_sep_presence(dBGraph* db_graph,
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

  Troolean result = is_sep_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("sep\t");
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
	  print_json_item("sep", "positive",  false  );
	}
      else
	{
	  print_json_item("sep", "negative",  false );
	}
      
    }
}



    
	Troolean is_seu_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/seu.fa");

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
	void print_seu_presence(dBGraph* db_graph,
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

  Troolean result = is_seu_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("seu\t");
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
	  print_json_item("seu", "positive",  false  );
	}
      else
	{
	  print_json_item("seu", "negative",  false );
	}
      
    }
}



    
	Troolean is_tsst1_positive(dBGraph* db_graph,
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
  strbuf_append_str(fa, "data/staph/virulence/tsst1.fa");

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
	void print_tsst1_presence(dBGraph* db_graph,
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

  Troolean result = is_tsst1_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("tsst1\t");
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
	  print_json_item("tsst1", "positive", true  );
	}
      else
	{
	  print_json_item("tsst1", "negative", true );
	}
      
    }
}




