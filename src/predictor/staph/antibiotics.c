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


void map_antibiotic_enum_to_str(Antibiotic ab, StrBuf* name)
{
  if (ab==NoDrug)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "NoDrug");
    }
  else if (ab==Gentamycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Gentamycin");
    }
  else if (ab==Penicillin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Penicillin");
    }
  else if (ab==Trimethoprim)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Trimethoprim");
    }
  else if (ab==Erythromycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Erythromycin");
    }
  else if (ab==Methicillin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Methicillin");
    }
  else if (ab==Ciprofloxacin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Ciprofloxacin");
    }
  else if (ab==Rifampicin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Rifampicin");
    }
  else if (ab==Tetracycline)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Tetracycline");
    }
  else if (ab==Mupirocin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Mupirocin");
    }
  else if (ab==FusidicAcid)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "FusidicAcid");
    }
  else if (ab==FusidicAcid)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "FusidicAcid");
    }
  else if (ab==Clindamycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Clindamycin");
    }
  else if (ab==Vancomycin)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "Vancomycin");
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
  abi->ab = Gentamycin;
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
  abi->ab = Penicillin;
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
  abi->ab = Trimethoprim;
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
  abi->ab = Erythromycin;
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
  abi->ab = Methicillin;
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
  abi->ab = Ciprofloxacin;
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
  abi->ab = Rifampicin;
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/rifampicin.fa");
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

boolean is_tetracycline_susceptible(dBGraph* db_graph,
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
  abi->ab = Tetracycline;
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/tetracycline.fa");
  abi->num_mutations = 0;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

 if (abi->genes[tetK]->percent_nonzero > GENE_THRESH_tetK)
    {
      return false;
    }  
 else if (abi->genes[tetL]->percent_nonzero > GENE_THRESH_tetL)
    {
      return false;
    }  
 else if (abi->genes[tetM]->percent_nonzero > GENE_THRESH_tetM)
    {
      return false;
    }  


 return true;
}


boolean is_mupirocin_susceptible(dBGraph* db_graph,
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
  abi->ab = Mupirocin;
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/mupirocin.fa");
  abi->num_mutations = 0;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

 if (abi->genes[mupA]->percent_nonzero > GENE_THRESH_mupA)
    {
      return false;
    }  
 else if (abi->genes[mupB]->percent_nonzero > GENE_THRESH_mupB)
    {
      return false;
    }  
 return true;
}


boolean is_fusidic_acid_susceptible(dBGraph* db_graph,
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
  abi->ab = FusidicAcid;
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/fusidic_acid.fa");
  abi->num_mutations = 915;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);


  if (abi->mut[fusA_A655E]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_B434N]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_E444K]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if ( 
	   (abi->mut[fusA_F652S]->some_resistant_allele_present==true)
	   &&
	   (abi->mut[fusA_Y654N]->some_resistant_allele_present==true)
	    )
    {
      return false;
    }
  else if (abi->mut[fusA_G451V]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_G452C]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_G452S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_G556S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_G617D]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_G664S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_H438N]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_H457Q]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_H457Y]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if ( 
	   (abi->mut[fusA_L461F]->some_resistant_allele_present==true)
	   &&
	   (abi->mut[fusA_A376V]->some_resistant_allele_present==true)
	   &&
	   (abi->mut[fusA_A655P]->some_resistant_allele_present==true)
	   &&
	   (abi->mut[fusA_D463G]->some_resistant_allele_present==true)
	    )
    {
      return false;
    }
  else if ( 
	   (abi->mut[fusA_L461F]->some_resistant_allele_present==true)
	   &&
	   (abi->mut[fusA_E444V]->some_resistant_allele_present==true)
	    )
    {
      return false;
    }
  else if (abi->mut[fusA_L461K]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_L461S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_M453I]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_M651I]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_P114H]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_P404L]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_P404Q]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_P406L]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_P478S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_Q115L]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R464C]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R464H]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R464S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R659C]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R659H]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R659L]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_R659S]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_T385N]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_T436I]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_T656K]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->mut[fusA_V90I]->some_resistant_allele_present==true)
    {
      return false;
    }
  else if (abi->genes[fusB]->percent_nonzero > GENE_THRESH_fusB)
    {
      return false;
    }  
  else if (abi->genes[fusC]->percent_nonzero > GENE_THRESH_fusC)
    {
      return false;
    }  

  else
    {
      return true;
    }
}


boolean is_clindamycin_susceptible(dBGraph* db_graph,
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
  //constitutuve only. inducible you get by checking erythromycin also,
  reset_antibiotic_info(abi);
  
  //setup antibiotic info object
  abi->ab = Clindamycin;
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/clindamycin.fa");
  abi->num_mutations = 0;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

 if (abi->genes[vga_A_LC]->percent_nonzero > GENE_THRESH_vga_A_LC)
    {
      return false;
    }  
 return true;

}


boolean is_vancomycin_susceptible(dBGraph* db_graph,
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
  //constitutuve only. inducible you get by checking erythromycin also,
  reset_antibiotic_info(abi);
  
  //setup antibiotic info object
  abi->ab = Vancomycin;
  strbuf_append_str(abi->fasta, "../data/staph/antibiotics/vancomycin.fa");
  abi->num_mutations = 0;

  load_antibiotic_mut_and_gene_info(db_graph,
				    file_reader,
				    abi,
				    rutils,
				    tmp_rvi,
				    tmp_gi);

 if (abi->genes[vanA]->percent_nonzero > GENE_THRESH_vanA)
    {
      return false;
    }  
 return true;

}



boolean print_antibiotic_susceptibility(dBGraph* db_graph,
					int (*file_reader)(FILE * fp, 
							   Sequence * seq, 
							   int max_read_length, 
							   boolean new_entry, 
							   boolean * full_entry),
					ReadingUtils* rutils,
					ResVarInfo* tmp_rvi,
					GeneInfo* tmp_gi,
					AntibioticInfo* abi,
					boolean (*func)(dBGraph* db_graph,
							int (*file_reader)(FILE * fp, 
									   Sequence * seq, 
									   int max_read_length, 
									   boolean new_entry, 
									   boolean * full_entry),
							ReadingUtils* rutils,
							ResVarInfo* tmp_rvi,
							GeneInfo* tmp_gi,
							AntibioticInfo* abi),
					StrBuf* tmpbuf
					)
{
  boolean suc;
  
  suc  = func(db_graph,
	      file_reader,
	      rutils,
	      tmp_rvi,
	      tmp_gi,
	      abi);

  map_antibiotic_enum_to_str(abi->ab, tmpbuf);
  printf("%s\t", tmpbuf->buff);
  if (suc==true)
    {
      printf("S\n");
    }
  else
    {
      printf("R\n");
    }
  return suc;
}


boolean print_clindamycin_susceptibility(dBGraph* db_graph,
					 int (*file_reader)(FILE * fp, 
							    Sequence * seq, 
							    int max_read_length, 
							    boolean new_entry, 
							    boolean * full_entry),
					 ReadingUtils* rutils,
					 ResVarInfo* tmp_rvi,
					 GeneInfo* tmp_gi,
					 AntibioticInfo* abi,
					 boolean (*func)(dBGraph* db_graph,
							 int (*file_reader)(FILE * fp, 
									    Sequence * seq, 
									    int max_read_length, 
									    boolean new_entry, 
									    boolean * full_entry),
							 ReadingUtils* rutils,
							 ResVarInfo* tmp_rvi,
							 GeneInfo* tmp_gi,
							 AntibioticInfo* abi),
					 StrBuf* tmpbuf,
					 boolean erythromycin_susceptible
					 )
{
  boolean suc;
  
  suc  = func(db_graph,
	      file_reader,
	      rutils,
	      tmp_rvi,
	      tmp_gi,
	      abi);

  map_antibiotic_enum_to_str(abi->ab, tmpbuf);
  printf("%s\t", tmpbuf->buff);
  if (suc==false)
    {
      printf("R(constitutive)\n");
    }
  else if (erythromycin_susceptible==false)
    {
      printf("R(inducible)\n");
    }
  else
    {
      printf("S\n");
    }
  return suc;
}

				  
