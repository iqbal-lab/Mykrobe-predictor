/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * **********************************************************************
 *
 * This file is part of myKrobe.
 *
 * myKrobe is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * myKrobe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with myKrobe.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
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
                  strbuf_append_str(name, "rifampicin");
                }
        
  else if (ab==isoniazid)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "isoniazid");
                }
        
  else if (ab==pyrazinamide)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "pyrazinamide");
                }
        
  else if (ab==ethambutol)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "ethambutol");
                }
        
  else if (ab==kanamycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "kanamycin");
                }
        
  else if (ab==capreomycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "capreomycin");
                }
        
  else if (ab==amikacin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "amikacin");
                }
        
  else if (ab==kanamycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "kanamycin");
                }
        
  else if (ab==streptomycin)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "streptomycin");
                }
        
  else if (ab==quinolones)
                {
                  strbuf_reset(name);
                  strbuf_append_str(name, "quinolones");
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
      abi->mut = (ResVarInfo**) malloc(sizeof(ResVarInfo*)*NUM_KNOWN_MUTATIONS);
      if (abi->mut==NULL)
	{
	  strbuf_free(abi->m_fasta);
	  free(abi);
	  return NULL;
	}
      abi->genes = (GeneInfo**) malloc(sizeof(GeneInfo*)*NUM_GENE_PRESENCE_GENES);
      if (abi->genes==NULL)
	{
	  free(abi->mut);
	  strbuf_free(abi->m_fasta);
	  free(abi);
	  return NULL;
	}
      abi->which_genes = (int*) calloc(MAX_GENES_PER_ANTIBIOTIC, sizeof(int));
      if (abi->which_genes==NULL)
	{
	  free(abi->genes);
	  free(abi->mut);
	  strbuf_free(abi->m_fasta);
	  free(abi);
	  return NULL;
	  
	}
      int i;
      for (i=0; i<NUM_KNOWN_MUTATIONS; i++)
	{
	  abi->mut[i] = alloc_and_init_res_var_info();
	}
      for (i=0; i<NUM_GENE_PRESENCE_GENES; i++)
	{
	  abi->genes[i] = alloc_and_init_gene_info();
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
	  free_res_var_info(abi->mut[i]);
	}
      for (i=0; i<NUM_GENE_PRESENCE_GENES; i++)
	{
	  free_gene_info(abi->genes[i]);
	}
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
      reset_res_var_info(abi->mut[i]);
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
					      ResVarInfo* tmp_rvi,	
					      int ignore_first, int ignore_last, int expected_covg)
{
  reset_reading_utils(rutils);
  reset_res_var_info(tmp_rvi);

  StrBuf* tmp1 = strbuf_new();
  StrBuf* tmp2 = strbuf_new();
  StrBuf* tmp3 = strbuf_new();

  
  int i;

  KnownMutation m = NotSpecified;
  for (i=0; i<abi->num_mutations; i++)
    {
      get_next_mutation_allele_info(fp, 
				    db_graph, 
				    tmp_rvi,
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

      copy_res_var_info(tmp_rvi, abi->mut[tmp_rvi->var_id]);
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
				       ResVarInfo* tmp_rvi,
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
					      tmp_rvi,
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






    Troolean is_amikacin_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = amikacin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/amikacin.fa");

      abi->num_mutations = 18;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rrs_A1401X;
      int last_mut = rrs_G1484X;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}


 
    Troolean is_streptomycin_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = streptomycin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/streptomycin.fa");

      abi->num_mutations = 47;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = gidB_A138V;
      int last_mut = rpsL_K88R;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}



    Troolean is_rifampicin_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = rifampicin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/rifampicin.fa");

      abi->num_mutations = 28;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rpoB_F425X;
      int last_mut = rpoB_L452X;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}




    Troolean is_quinolones_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = quinolones;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/quinolones.fa");

      abi->num_mutations = 60;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = gyrA_H85X;
      int last_mut = gyrA_D94X;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}



    Troolean is_pyrazinamide_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = pyrazinamide;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/pyrazinamide.fa");

      abi->num_mutations = 50;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = pncA_D49N;
      int last_mut = pncA_V21G;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}



    Troolean is_kanamycin_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = kanamycin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/kanamycin.fa");

      abi->num_mutations = 23;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rrs_A1401X;
      int last_mut = eis_Cu10T;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}


    Troolean is_isoniazid_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
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
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = katG_S315X;
      int last_mut = fabG1_Au16X;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}



    Troolean is_ethambutol_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = ethambutol;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/ethambutol.fa");

      abi->num_mutations = 18;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = embB_M306X;
      int last_mut = embB_G406S;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}




    Troolean is_capreomycin_susceptible(dBGraph* db_graph,
               int (*file_reader)(FILE * fp, 
                          Sequence * seq, 
                          int max_read_length, 
                          boolean new_entry, 
                          boolean * full_entry),
                 ReadingUtils* rutils,
                 ResVarInfo* tmp_rvi,
                 GeneInfo* tmp_gi,
                 AntibioticInfo* abi,
                 StrBuf* install_dir,
                 int ignore_first, int ignore_last, int expected_covg,
                 double lambda_g, double lambda_e, double err_rate)
    {
      reset_antibiotic_info(abi);
      
      //setup antibiotic info object
      abi->ab = capreomycin;
      strbuf_append_str(abi->m_fasta, install_dir->buff);

      strbuf_append_str(abi->m_fasta, "data/tb/antibiotics/capreomycin.fa");

      abi->num_mutations = 18;
      abi->num_genes=0;
      double epsilon = pow(1-err_rate, db_graph->kmer_size);
      load_antibiotic_mut_and_gene_info(db_graph,
                        file_reader,
                        abi,
                        rutils,
                        tmp_rvi,
                        tmp_gi,
                        ignore_first, ignore_last, expected_covg,
        install_dir);




      int first_mut = rrs_A1401X;
      int last_mut = rrs_G1484X;

      int i;
       //if you have any of these resistance alleles - call resistant
  double min_conf=9999999999;
  for (i=first_mut; i<=last_mut; i++)
    {
      Model best_model;
      InfectionType I=
    resistotype(abi->mut[i],
           err_rate, db_graph->kmer_size, lambda_g, lambda_e, epsilon,
            &best_model, MaxLikelihood);
      if (min_conf>best_model.conf)
    {
      min_conf=best_model.conf;
    }
      if (I==Resistant)
    {
      return _False;
    }
    }


  if (min_conf>MIN_CONFIDENCE)
    {
      return _True;
    }
  else
    {
      return _Inconclusive;
    }

}


void print_antibiotic_susceptibility(dBGraph* db_graph,
					int (*file_reader)(FILE * fp, 
							   Sequence * seq, 
							   int max_read_length, 
							   boolean new_entry, 
							   boolean * full_entry),
					ReadingUtils* rutils,
					ResVarInfo* tmp_rvi,
					GeneInfo* tmp_gi,
					AntibioticInfo* abi,
					Troolean (*func)(dBGraph* db_graph,
							int (*file_reader)(FILE * fp, 
									   Sequence * seq, 
									   int max_read_length, 
									   boolean new_entry, 
									   boolean * full_entry),
							ReadingUtils* rutils,
							ResVarInfo* tmp_rvi,
							GeneInfo* tmp_gi,
							AntibioticInfo* abi,
							StrBuf* install_dir,
							int ignore_first, int ignore_last, int expected_covg,
							double lambda_g, double lambda_e, double err_rate),
					StrBuf* tmpbuf,
					StrBuf* install_dir,
					int ignore_first, int ignore_last,
					int expected_covg,
					double lambda_g, double lambda_e, double err_rate,
					OutputFormat format,
					boolean output_last//for JSON
					)
{
  Troolean suc;
  
  suc  = func(db_graph,
	      file_reader,
	      rutils,
	      tmp_rvi,
	      tmp_gi,
	      abi, 
	      install_dir,
	      ignore_first, 
	      ignore_last, 
	      expected_covg,
	      lambda_g,
	      lambda_e,
	      err_rate);

  
  map_antibiotic_enum_to_str(abi->ab, tmpbuf);
  if (format==Stdout)
    {
      printf("%s\t", tmpbuf->buff);
      if (suc==_True)
	{
	  printf("S\n");
	}
      else if (suc==_False)
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
      if (suc==_True)
	{
	    print_json_item(tmpbuf->buff, "S", output_last);
	}
      else if (suc==_False)
	{
	  print_json_item(tmpbuf->buff, "R", output_last);
	}
      else
	{
	  print_json_item(tmpbuf->buff, "Inconclusive", output_last);
	}
    }

}
