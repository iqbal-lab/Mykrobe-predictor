#include "global.h" // Covg def 
#include "seq.h" // Sequence def 
#include "binary_kmer.h" // KmerSlidingWindow
#include "dB_graph.h"// dBGraph
#include "db_complex_genotyping.h" // CovgArray

typedef struct
{
  Covg median_covg;
  Covg min_covg;
  int  percent_nonzero;
  StrBuf* strbuf;
  StrBuf* name;
} GeneInfo;

GeneInfo* alloc_and_init_gene_info();
void free_gene_info(GeneInfo* rvi);
void reset_gene_info(GeneInfo* rvi);

void get_next_gene_info(FILE* fp, 
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
