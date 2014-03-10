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


AntibioticInfo* alloc_antibiotic_info()
{
  AntibioticInfo* abi = (AntibioticInfo*) calloc(1, sizeof(AntibioticInfo));
  if (abi==NULL)
    {
      return NULL;
    }
  else
    {
      abi->fasta = strbuf_new();
      int i;
      abi->mut = (ResVarInfo**) malloc(sizeof(ResVarInfo*)*NUM_KNOWN_MUTATIONS);
      if (abi->mut==NULL)
	{
	  strbuf_free(abi->fasta);
	  free(abi);
	  return NULL;
	}
      abi->genes = (GeneInfo**) malloc(sizeof(GeneInfo*)*NUM_GENE_PRESENCE_GENES);
      if (abi->genes==NULL)
	{
	  free(abi->mut);
	  strbuf_free(abi->fasta);
	  free(abi);
	  return NULL;
	}
      
      for (i=0; i<NUM_KNOWN_MUTATIONS; i++)
	{
	  abi->mut[i] = alloc_and_init_res_var_info();
	}
      for (i=0; i<NUM_GENE_PRESENCE_GENES; i++)
	{
	  abi->genes[i] = alloc_and_init_gene_info();
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
      strbuf_free(abi->fasta);
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
  strbuf_reset(abi->fasta);
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
					      ResVarInfo* tmp_rvi)
{
  reset_reading_utils(rutils);
  reset_res_var_info(tmp_rvi);

  StrBuf* tmp1 = strbuf_new();
  StrBuf* tmp2 = strbuf_new();
  StrBuf* tmp3 = strbuf_new();

  
  int i;

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
				    tmp1, tmp2, tmp3);

      copy_res_var_info(tmp_rvi, abi->mut[tmp_rvi->var_id]);
    }
  
  strbuf_free(tmp1);
  strbuf_free(tmp2);
  strbuf_free(tmp3);

}



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

      copy_gene_info(tmp_gi, abi->genes[tmp_gi->name]);
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
				       GeneInfo* tmp_gi)
{

  FILE* fp = fopen(abi->fasta->buff, "r");
  if (fp==NULL)
    {
      die("Cannot open %s - should be there as part of the install - did you run out of disk mid-install?\n",
	  abi->fasta->buff);
    }
  
  load_antibiotic_mutation_info_on_sample(fp,
					  db_graph,
					  file_reader,
					  abi,
					  rutils, 
					  tmp_rvi);
  load_antibiotic_gene_presence_info_on_sample(fp,
					       db_graph,
					       file_reader,
					       abi,
					       rutils,
					       tmp_gi);
  fclose(fp);
}






boolean is_gentamycin_susceptible(dBGraph* db_graph,
				  int (*file_reader)(FILE * fp, 
						     Sequence * seq, 
						     int max_read_length, 
						     boolean new_entry, 
						     boolean * full_entry),
				  ReadingUtils* rutils,
				  ResVarInfo* tmp_rvi,
				  GeneInfo* tmp_gi,
				  AntibioticInfo* abi
				  )
{
  reset_antibiotic_info(abi);

  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/gentamycin.fa");
  abi->num_mutations = 0;//entirely determined by gene presence

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

  if (abi->genes[aacAaphD]->percent_nonzero > GENE_THRESH_aacAaphD)
    {
      return false;
    }
  return true;
  
}


boolean is_penicillin_susceptible(dBGraph* db_graph,
				  int (*file_reader)(FILE * fp, 
						     Sequence * seq, 
						     int max_read_length, 
						     boolean new_entry, 
						     boolean * full_entry),
				  ReadingUtils* rutils,
				  ResVarInfo* tmp_rvi,
				  GeneInfo* tmp_gi,
				  AntibioticInfo* abi
				  )

{

  reset_antibiotic_info(abi);

  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/penicillin.fa");
  abi->num_mutations = 0;//entirely determined by gene presence

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

  if (abi->genes[blaZ]->percent_nonzero > GENE_THRESH_blaZ)
    {
      return false;
    }
  return true;

}


boolean is_trimethoprim_susceptible(dBGraph* db_graph,
				    int (*file_reader)(FILE * fp, 
						       Sequence * seq, 
						       int max_read_length, 
						       boolean new_entry, 
						       boolean * full_entry),
				    ReadingUtils* rutils,
				    ResVarInfo* tmp_rvi,
				    GeneInfo* tmp_gi,
				    AntibioticInfo* abi
				  )

{

  reset_antibiotic_info(abi);

  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/trimethoprim.fa");
  abi->num_mutations = 113;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

  if (abi->mut[dfrB_F99I]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[dfrB_F99S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[dfrB_F99Y]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[dfrB_H31N]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[dfrB_L41F]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->genes[dfrG]->percent_nonzero > GENE_THRESH_blaZ)
    {
      return false;
    }
  else
    {
      return true;
    }
}

boolean is_erythromycin_susceptible(dBGraph* db_graph,
				    int (*file_reader)(FILE * fp, 
						       Sequence * seq, 
						       int max_read_length, 
						       boolean new_entry, 
						       boolean * full_entry),
				    ReadingUtils* rutils,
				    ResVarInfo* tmp_rvi,
				    GeneInfo* tmp_gi,
				    AntibioticInfo* abi)
{
  reset_antibiotic_info(abi);
  
  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/erythromycin.fa");
  abi->num_mutations = 0;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

 if (abi->genes[ermA]->percent_nonzero > GENE_THRESH_ermA)
    {
      return false;
    }
 else if (abi->genes[ermB]->percent_nonzero > GENE_THRESH_ermB)
    {
      return false;
    }
 else if (abi->genes[ermC]->percent_nonzero > GENE_THRESH_ermC)
    {
      return false;
    }
 else if (abi->genes[ermT]->percent_nonzero > GENE_THRESH_ermT)
    {
      return false;
    }
 else if (abi->genes[msrA]->percent_nonzero > GENE_THRESH_msrA)
    {
      return false;
    }
 else
   {
     return true;
   }
}


boolean is_methicillin_susceptible(dBGraph* db_graph,
				   int (*file_reader)(FILE * fp, 
						      Sequence * seq, 
						      int max_read_length, 
						      boolean new_entry, 
						      boolean * full_entry),
				   ReadingUtils* rutils,
				   ResVarInfo* tmp_rvi,
				   GeneInfo* tmp_gi,
				   AntibioticInfo* abi)
{
  reset_antibiotic_info(abi);
  
  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/methicillin.fa");
  abi->num_mutations = 0;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

 if (abi->genes[mecA]->percent_nonzero > GENE_THRESH_mecA)
    {
      return false;
    }  
 return true;
}


boolean is_ciprofloxacin_susceptible(dBGraph* db_graph,
				   int (*file_reader)(FILE * fp, 
						      Sequence * seq, 
						      int max_read_length, 
						      boolean new_entry, 
						      boolean * full_entry),
				   ReadingUtils* rutils,
				   ResVarInfo* tmp_rvi,
				   GeneInfo* tmp_gi,
				   AntibioticInfo* abi)
{
  reset_antibiotic_info(abi);
  
  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/ciprofloxacin.fa");
  abi->num_mutations = 94;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

  if (abi->mut[gyrA_E88K]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[gyrA_S84A]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[gyrA_S84L]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[gyrA_S85P]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[grlA_S80F]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[grlA_S80Y]->some_resistant_allele_present==true)
    {
      return false;
    }
  else
    {
      return true;
    }
}


boolean is_rifampicin_susceptible(dBGraph* db_graph,
				   int (*file_reader)(FILE * fp, 
						      Sequence * seq, 
						      int max_read_length, 
						      boolean new_entry, 
						      boolean * full_entry),
				   ReadingUtils* rutils,
				   ResVarInfo* tmp_rvi,
				   GeneInfo* tmp_gi,
				   AntibioticInfo* abi)
{
  reset_antibiotic_info(abi);
  
  //setup antibiotic info object
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/ciprofloxacin.fa");
  abi->num_mutations = 430;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);
  if (abi->mut[rpoB_A477D]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[rpoB_A477V]->some_resistant_allele_present==true)
    {
      return false;
    } 



  else if (abi->mut[rpoB_D471Y]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_D550G]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_H481D]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_H481N]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_H481Y]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_I527F]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_ins475G]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_ins475H]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if ( 
	   (abi->mut[rpoB_M470T]->some_resistant_allele_present==true)
	   &&
	   (abi->mut[rpoB_D471G]->some_resistant_allele_present==true)
	    )
    {
      return false;
    } 
  else if (abi->mut[rpoB_Q456K]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_Q468K]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_Q468L]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_Q468R]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_R484H]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_S463P]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_S464P]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else if (abi->mut[rpoB_S486L]->some_resistant_allele_present==true)
    {
      return false;
    } 
  else
    {
      return true;
    }

}
/*
boolean is_tetracycline_susceptible(dBGraph* db_graph)
{
}

boolean is_mupirocin_susceptible(dBGraph* db_graph)
{
}

boolean is_fusidic_acid_susceptible(dBGraph* db_graph)
{
}

boolean is_clindamycin_susceptible(dBGraph* db_graph)
{
//constitutuve only. inducible you get by checking erythromycin also,
}

*/
				  
