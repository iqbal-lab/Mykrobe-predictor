/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 *  json.c 
*/


#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "myk_pred_global.h"

void print_json_start()
{
  fflush(stdout);
  printf("{\n");
}

void print_json_version()
{
  printf("\t\"version\":  \"%d.%d.%d\",\n", MYK_VERSION,MYK_SUBVERSION,MYK_SUBSUBVERSION);
}

void print_json_end()
{
  printf("}\n");
  fflush(stdout);
}

// void print_json_species_start()
// {
//   printf("\t\t\"species\": {\n");
// }
// void print_json_phylo_group_start()
// {
//   printf("\t\t\"phylo_group\": {\n");
// }

// void print_json_lineage_start()
// {
//   printf("\t\t\"lineage\": {\n");
// }


// void print_json_species_end()
// {
//   printf("\t\t},\n");
// }
void print_json_phylo_group_end()
{
  printf("\t\t},\n");
}
// void print_json_lineage_end()
// {
//   printf("\t\t}\n");
// }

void print_json_phylogenetics_start()
{
  printf("\t\"phylogenetics\": {\n");
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

void print_json_called_variants_end(boolean last)
{
  if (last==true)
  {
    printf("\t}\n");
  }
  else{
    printf("\t},\n");
  }
  
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

void print_json_called_end(boolean last)
{
  if (last==false){
    printf("\t\t},\n");
  }
  else{
    printf("\t\t}\n");
  }
  
}

void print_json_called_variant_end(boolean last)
{
  print_json_called_end(last);
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

void print_json_called_gene_end(boolean last)
{
  print_json_called_end(last);
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


