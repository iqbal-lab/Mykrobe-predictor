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
  gene_presence.h - build the de Bruijn graph
*/

#ifndef GENE_PRESENCE_H
#define GENE_PRESENCE_H


#include "global.h" // Covg def 
#include "seq.h" // Sequence def 
#include "binary_kmer.h" // KmerSlidingWindow
#include "dB_graph.h"// dBGraph
#include "dB_graph_supernode.h"




typedef enum
  {
    aacAaphD=0,
    blaZ=1,
    dfrA=2,
    dfrG=3,
    ermA=4,
    ermB=5,
    ermC=6,
    ermT=7,
    fusB=8,
    fusC=9,
    vga_A_LC=10,
    msrA=11,
    mecA=12,
    tetK=13,
    tetL=14,
    tetM=15,
    vanA=16,
    mupA=17,
    mupB=18,
    luk=19,//lukS-PV and lukF-PV
    unspecified_gpg = 20,
  } GenePresenceGene;

#define NUM_GENE_PRESENCE_GENES 20    //ignore unspecified_gpg
#define MAX_LEN_GENE 3110

GenePresenceGene map_string_to_gene_presence_gene(StrBuf* sbuf);
boolean map_gene_to_fasta(GenePresenceGene gene, StrBuf* fa, StrBuf* install_dir);

typedef struct
{
  Covg median_covg;
  Covg min_covg;
  Covg median_covg_on_nonzero_nodes;
  int  percent_nonzero;
  int num_instances;//will read several exemplars of this gene and store info on the best match
  StrBuf* strbuf;
  GenePresenceGene name;
} GeneInfo;

GeneInfo* alloc_and_init_gene_info();
void free_gene_info(GeneInfo* rvi);
void reset_gene_info(GeneInfo* rvi);
void copy_gene_info(GeneInfo* from_gi, GeneInfo* to_gi);

int get_next_gene_info(FILE* fp, 
		       dBGraph* db_graph, 
		       GeneInfo* gene_info,
		       Sequence* seq, 
		       KmerSlidingWindow* kmer_window,
		       int (*file_reader)(FILE * fp, 
					  Sequence * seq, 
					  int max_read_length, 
					  boolean new_entry, 
					  boolean * full_entry),
		       dBNode** array_nodes, 
		       Orientation*  array_or,
		       CovgArray* working_ca, int max_read_length);
const char* map_enum_to_gene_name(GenePresenceGene gene);


#endif
