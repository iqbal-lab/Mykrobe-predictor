/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  json.h
*/

#ifndef JSON_H_
#define JSON_H_


typedef enum
  {
    Stdout = 0,
    JSON = 1,
  } OutputFormat;

void print_json_start();
void print_json_end();

// void print_json_species_start();
// void print_json_species_end();
// void print_json_phylo_group_start();
void print_json_phylo_group_end();
// void print_json_lineage_start();
// void print_json_lineage_end();


void print_json_phylogenetics_start();
void print_json_phylogenetics_end();

void print_json_susceptibility_start();
void print_json_susceptibility_end();
void print_json_called_variants_start();
void print_json_called_variants_end();
void print_json_called_variant_start(const char* str1);
void print_json_called_variant_item(char* str1, int val, boolean last);
void print_json_called_variant_end(boolean last);

void print_json_called_genes_start();
void print_json_called_genes_end();
void print_json_called_gene_start(StrBuf* sbuf);
void print_json_called_gene_item(char* str1, int val, boolean last);
void print_json_called_gene_end(boolean last);
void print_json_called_end(boolean last);
void print_json_virulence_start();
void print_json_virulence_end();
void print_json_item(char* str1, char* str2, boolean output_last);

#endif

