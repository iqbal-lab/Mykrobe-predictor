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
  printf("\tSpecies:{\n");
}

void print_json_species_end()
{
  printf("\t},\n");
}

void print_json_susceptibility_start()
{
  printf("\tSusceptibility:{\n");
}

void print_json_susceptibility_end()
{
  printf("\t},\n");
}

void print_json_next_item(char* str)
{
  printf("%s,\n", str);
}

void print_json_last_item(char* str)
{
  printf("%s\n", str);
}



