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
  json.h
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
void print_json_species_start();
void print_json_species_end();
void print_json_lineage_start();
void print_json_lineage_end();
void print_json_phylogenetics_start();
void print_json_phylogenetics_end();

void print_json_susceptibility_start();
void print_json_susceptibility_end();
void print_json_called_variants_start();
void print_json_called_variants_end();
void print_json_called_variant_start(const char* str1);
void print_json_called_variant_item(char* str1, int val, boolean last);
void print_json_called_variant_end();

void print_json_called_genes_start();
void print_json_called_genes_end();
void print_json_called_gene_start(const char* str1);
void print_json_called_gene_item(char* str1, int val, boolean last);
void print_json_called_gene_end();

void print_json_virulence_start();
void print_json_virulence_end();
void print_json_item(char* str1, char* str2, boolean output_last);

#endif

