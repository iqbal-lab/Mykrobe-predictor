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
  antibiotics.h
*/

#ifndef ANTIBIOTICS_H_
#define ANTIBIOTICS_H_

#include "stdlib.h"
#include "string_buffer.h"
#include "gene_presence.h"
#include "genotyping_known.h"
#include "known_mutations.h"
#include "file_reader.h"



//#define MAX_MUTS_IN_ANY_GENE 129
#define MAX_LEN_MUT_ALLELE 61


#define GENE_THRESH_blaZ 25
#define GENE_THRESH_msrA 90
#define GENE_THRESH_aacAaphD 70
#define GENE_THRESH_ermC 80
#define GENE_THRESH_tetK 80
#define GENE_THRESH_dfrA 70 
#define GENE_THRESH_dfrG 70
#define GENE_THRESH_ermA 70
#define GENE_THRESH_ermB 70
#define GENE_THRESH_ermT 70
#define GENE_THRESH_fusB 70
#define GENE_THRESH_fusC 70
#define GENE_THRESH_vga_A_LC 70



typedef enum 
 {
   NoDrug=0,
   Gentamycin=1,
   Penicillin=2,
   Trimethoprim=3,
   Erythromycin=4,
   Methicillin=5,
   Ciprofloxacin=6,
   Rifampicin=7,
   Tetracycline=8,
   Mupirocin=9,
   FusidicAcid=10,
   Clindamycin=11,
  } Antibiotic;





typedef struct
{
  Antibiotic ab;
  StrBuf* fasta; //mutations and genes  
  int num_mutations;
  ResVarInfo** mut;// array of length NUM_KNOWN_MUTATIONS];//indexed by enum of mutation names
  GeneInfo** genes;// array of length NUM_GENE_PRESENCE_GENES];//indexed by enum of gene names
  //which genes to check for presence of is only held implicitly, by reading
  //the relevant file.
  
} AntibioticInfo;






AntibioticInfo* alloc_antibiotic_info();
void free_antibiotic_info(AntibioticInfo* abi);
void reset_antibiotic_info(AntibioticInfo* abi);

void  load_antibiotic_mutation_info_on_sample(FILE* fp,
					      dBGraph* db_graph,
					      int (*file_reader)(FILE * fp, 
								 Sequence * seq, 
								 int max_read_length, 
								 boolean new_entry, 
								 boolean * full_entry),
					      AntibioticInfo* abi,
					      ReadingUtils* rutils,
					      ResVarInfo* tmp_rvi);

void load_antibiotic_gene_presence_info_on_sample(FILE* fp,
						  dBGraph* db_graph,
						  int (*file_reader)(FILE * fp, 
								     Sequence * seq, 
								     int max_read_length, 
								     boolean new_entry, 
								     boolean * full_entry),
						  AntibioticInfo* abi,
						  ReadingUtils* rutils,
						  GeneInfo* tmp_gi);

void load_antibiotic_mut_and_gene_info(dBGraph* db_graph,
				       int (*file_reader)(FILE * fp, 
							  Sequence * seq, 
							  int max_read_length, 
							  boolean new_entry, 
							  boolean * full_entry),
				       AntibioticInfo* abi,
				       ReadingUtils* rutils,
				       ResVarInfo* tmp_rvi,
				       GeneInfo* tmp_gi);

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
				  );


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
				  );

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
				  );


#endif
