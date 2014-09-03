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
  json.c 
*/


#include <stdio.h>
#include <stdlib.h>
#include "global.h"

void print_json_start()
{
  printf("{\n");
}

void print_json_end()
{
  printf("}\n");
}

void print_json_species_start()
{
  printf("\t\t\"species\": {\n");
}
void print_json_phylogenetics_start()
{
  printf("\t\"phylogenetics\": {\n");
}
void print_json_lineage_start()
{
  printf("\t\t\"lineage\": {\n");
}


void print_json_species_end()
{
  printf("\t\t},\n");
}
void print_json_lineage_end()
{
  printf("\t\t},\n");
}
void print_json_phylogenetics_end()
{
  printf("\t},\n");
}

void print_json_susceptibility_start()
{
  printf("\t\"susceptibility\" :{\n");
}

void print_json_susceptibility_end()
{
  printf("\t},\n");
}

// Json printing of variants
void print_json_called_variants_start()
{
  printf("\t\"called_variants\" :{\n");
}
void print_json_called_variants_end()
{
  printf("\t},\n");
}


void print_json_called_variant_start(const char* str1)
{
  printf("\t\t\"%s\" :{\n",str1);
}
void print_json_called_variant_item(char* str1, int val, boolean last)
{
  printf("\t\t\t\"%s\": \"%i\"", str1, val);
  if (last==false)
    {
      printf(",");
    }
  printf("\n");
}

void print_json_called_variant_end()
{
  printf("\t\t},\n");
}
//
// Json printing of genes
void print_json_called_genes_start()
{
  printf("\t\"called_genes\" :{\n");
}
void print_json_called_genes_end()
{
  printf("\t},\n");
}


void print_json_called_gene_start(const char* str1)
{
  printf("\t\t\"%s\" :{\n",str1);
}
void print_json_called_gene_item(char* str1, int val, boolean last)
{
  printf("\t\t\t\"%s\": \"%i\"", str1, val);
  if (last==false)
    {
      printf(",");
    }
  printf("\n");
}

void print_json_called_gene_end()
{
  printf("\t\t},\n");
}
//

void print_json_virulence_start()
{
  printf("\t\"virulence_toxins\" :{\n");
}

void print_json_virulence_end()
{
  printf("\t}\n");
}


void print_json_item(char* str1, char* str2, boolean last)
{
  printf("\t\t\"%s\": \"%s\"", str1, str2);
  if (last==false)
    {
      printf(",");
    }
  printf("\n");
}


