/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  file_reader.h
*/

#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <sys/stat.h>

#include <string_buffer.h>

#include "global.h"
#include "seq.h"
#include "dB_graph.h"
#include "file_format.h"
#include "graph_info.h"
#include "db_variants.h"

extern int MAX_FILENAME_LENGTH;
extern int MAX_READ_LENGTH;

typedef enum {
  EValid                        = 0,
  ECannotReadMagicNumber        = 1,
  ECanReadMagicNumberButIsWrong = 2,
  ECannotReadBinversion         = 3,
  EValidButOldBinVersion        = 4,
  EInvalidBinversion            = 5,
  ECannotReadKmer               = 6,
  EWrongKmer                    = 7,
  ECannotReadNumBitfields       = 8,
  EWrongNumberBitfields         = 9,
  ECannotReadNumColours         = 10,
  EBadColours                   = 11,
  EFailedToReadReadLensAndCovgs = 12,
  ECannotReadEndOfHeaderMagicNumber = 13,
  EFailedToReadSampleIds        =14,
  EFailedToReadSampleIdsSeemsTooLong = 15,
  EFailedToReadSeqErrRates      =16,
  EFailedToReadErrorCleaningInfo=17,
  ECanReadEndOfHeaderMagicNumberButIsWrong = 18,
  EGarbage = 19,
  EBinaryHasTooManyColoursGivenFirstColour=20, //ie starting at colour 10, loading 100 colours but 110>NUMBER_OF_COLOURS-1
} BinaryHeaderErrorCode;
typedef struct {
  int version;
  int kmer_size;
  int number_of_bitfields;
  int number_of_colours;
  GraphInfo* ginfo;
} BinaryHeaderInfo;

boolean dir_exists(char* dir_to_check);

// mkpath - ensure all directories in path exist
// Returns 1 on success, 0 on failure
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char mkpath(const char *path, mode_t mode);

StrBuf* file_reader_get_strbuf_of_dir_path(char* path);

boolean subsample_null();

void load_se_seq_data_into_graph_colour(
  const char *file_path,
  char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_se,
  char ascii_fq_offset, int colour_index, dBGraph *db_graph,
  unsigned long long *bad_reads, unsigned long long *dup_reads,
  unsigned long long *bases_read, unsigned long long *bases_loaded,
  unsigned long *readlen_count_array, unsigned long readlen_count_array_size,
  boolean (*subsample_func)() );

void load_pe_seq_data_into_graph_colour(
  const char *file_path1, const char *file_path2,
  char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_pe,
  char ascii_fq_offset, int colour_index, dBGraph *db_graph,
  unsigned long long *bad_reads, unsigned long long *dup_reads,
  unsigned long long *bases_read, unsigned long long *bases_loaded,
  unsigned long *readlen_count_array, unsigned long readlen_count_array_size,
  boolean (*subsample_func)() );

void load_se_filelist_into_graph_colour(
  char* se_filelist_path,
  int qual_thresh, int homopol_limit, boolean remove_dups_se,
  char ascii_fq_offset, int colour, dBGraph* db_graph, char is_colour_list,
  unsigned int *total_files_loaded,
  unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
  unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
  unsigned long *readlen_count_array, unsigned long readlen_count_array_size,
  boolean (*subsample_func)() );

void load_pe_filelists_into_graph_colour(
  char* pe_filelist_path1, char* pe_filelist_path2,
  int qual_thresh, int homopol_limit, boolean remove_dups_pe,
  char ascii_fq_offset, int colour, dBGraph* db_graph, char is_colour_lists,
  unsigned int *total_file_pairs_loaded,
  unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
  unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
  unsigned long *readlen_count_array, unsigned long readlen_count_array_size,
  boolean (*subsample_func)() );

// End of loading sequence data

void initialise_binary_header_info(BinaryHeaderInfo* binfo, GraphInfo* ginfo);


void  load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(
  KmerSlidingWindowSet * windows, boolean* prev_full_ent, //boolean* full_ent,
  long long* bases_loaded, boolean mark_read_starts, dBGraph* db_graph,
  int index, long long** read_len_count_array);


//pass in a single kmer sliding window and the Sequence* it was derived from. Will find the nodes correspinding to this seqeunce
//and put them in array. Also will check that edges exist as expected from the Sequence*
void load_kmers_from_sliding_window_into_array(KmerSlidingWindow* kmer_window, dBGraph* db_graph, dBNode** array_nodes, Orientation* array_orientations, 
					       int max_array_size, 
					       boolean require_nodes_to_lie_in_given_colour, int colour);


//use preallocated sliding window, and get all the kmers from the passed-in sequence. Any kmer that would have contained an N is returned as NULL
int get_single_kmer_sliding_window_from_sequence(char * seq, int length, short kmer_size, KmerSlidingWindow* kmer_window, dBGraph* db_graph);





// gets the next number_of_bases_to_load bases from fasta file, and returns them in the array of nodes.
// assumes this file has already been loaded into the graph.
// returns the number of nodes loaded. If this is less than what you asked for, you know it has hit the end of the file.
// We expect this to be used as follows:
// repeated calls of this function load etc into the LAST number_of_bases_to_load places of the relevant arrays

int load_seq_into_array(FILE* chrom_fptr, int number_of_nodes_to_load, int length_of_arrays, 
			dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry, dBGraph * db_graph);





//functions for loading multicolour graphs
long long load_multicolour_binary_from_filename_into_graph(char* filename,  dBGraph* db_graph, GraphInfo* ginfo, int* num_cols_in_loaded_binary); 

// last arg is to load the "union" of graphs. This is not a special case of load_multicolour
long long load_single_colour_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, 
								  GraphInfo* ginfo, 
								  boolean all_entries_are_unique, 
								  int colour_loading_into,
								  boolean only_load_kmers_already_in_hash, int colour_clean,
								  boolean load_all_kmers_but_only_increment_covg_on_new_ones);
								  

long long load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(char* filename,  dBGraph* db_graph, GraphInfo* db_graph_info,
											   boolean all_entries_are_unique, int index,
											   boolean only_load_kmers_already_in_hash, int colour_clean,
											   boolean load_all_kmers_but_only_increment_covg_on_new_ones);

long long load_population_as_binaries_from_graph(char* filename, int first_colour,boolean about_to_load_first_binary_into_empty_graph, 
						 dBGraph* db_graph, GraphInfo* db_graph_info,
						 boolean only_load_kmers_already_in_hash, int colour_clean,
						 boolean load_all_kmers_but_only_increment_covg_on_new_ones);

void dump_successive_cleaned_binaries(char* filename, int in_colour, int clean_colour, char* suffix, dBGraph* db_graph, GraphInfo* db_graph_info );



//functions for comparing graph with reference, or comparing reads with the graph

void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                  long long * bad_reads, int max_read_length, dBGraph * db_graph,
                                                  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);


/*
// DEV: flagged for removal -- not used, uses old file reader
void read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(
  FILE* fp, FILE* fout,
  int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,
                      boolean new_entry, boolean * full_entry), 
  long long * bad_reads, int max_read_length, dBGraph * db_graph, 
  int index, boolean is_for_testing,
  char** for_test_array_of_clean_reads, int* for_test_index);
*/

int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph);

int get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										     KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph, int index); 

void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);




//gets number_of_bases_to_load's worth of kmers, and returns the corresponding nodes, orientations etc in he array passed in.

int load_seq_into_array(FILE* chrom_fptr, int number_of_bases_to_load, int length_of_arrays,
			    dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			    Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry,  dBGraph * db_graph);


int align_next_read_to_graph_and_return_node_array(FILE* fp, int max_read_length, dBNode** array_nodes, Orientation* array_orientations, 
						   boolean require_nodes_to_lie_in_given_colour, boolean* full_entry,
						   int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
						   Sequence* seq, KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour);


int given_prev_kmer_align_next_read_to_graph_and_return_node_array_including_overlap(char* prev_kmer, FILE* fp, int max_read_length, 
										     dBNode** array_nodes, Orientation* array_orientations, 
										     boolean require_nodes_to_lie_in_given_colour,
										     boolean* full_entry,
										     int (* file_reader)(FILE * fp, Sequence * seq, 
													 int max_read_length,boolean new_entry, 
													 boolean * full_entry), 
										     Sequence* seq, Sequence* seq_inc_prev_kmer, 
										     KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour);

int read_next_variant_from_full_flank_file(FILE* fptr, int max_read_length,
					   VariantBranchesAndFlanks* var, dBGraph* db_graph, 
					   int (file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
					   Sequence* seq, Sequence* seq_inc_prev_kmer, KmerSlidingWindow* kmer_window);


//void print_binary_signature(FILE * fp,int kmer_size, int num_cols, int* array_mean_readlens, long long* array_total_seq);
void print_binary_signature_NEW(FILE * fp,int kmer_size, int num_cols, GraphInfo* ginfo, int first_col, int version);

//boolean check_binary_signature(FILE * fp,int kmer_size, int bin_version, int* number_of_colours_in_binary, int** array_mean_readlens, long long** array_total_seqs, int *return_binversion);
boolean check_binary_signature_NEW(FILE * fp,int kmer_size, 
				   BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode,
				   int first_colour_loading_into);

boolean query_binary_NEW(FILE * fp, BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode,
			 int first_colour_loading_into);

void print_error_cleaning_object(FILE* fp, GraphInfo* ginfo, int colour);

boolean get_extra_data_from_header(FILE * fp, BinaryHeaderInfo* binfo, 
				   BinaryHeaderErrorCode* ecode, int first_colour_loading_into);
boolean get_read_lengths_and_total_seqs_from_header(FILE * fp, BinaryHeaderInfo* binfo, 
						    BinaryHeaderErrorCode* ecode, int first_colour_loading_into);

boolean  get_binversion6_extra_data(FILE * fp, BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode, int first_colour_loading_into);
boolean read_next_error_cleaning_object(FILE* fp, ErrorCleaning* cl);


int load_paths_from_filelist(char* filelist_path, char** path_array);

boolean check_colour_list(char* filename, int kmer);
boolean check_ctx_list(char* filename, int kmer);

#endif /* FILE_READER_H_ */
