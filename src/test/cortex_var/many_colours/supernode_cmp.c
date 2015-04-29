/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 */
/*
  supernode_cmp.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int supernode_cmp(const void* supernode1, const void* supernode2)
{
  return strcmp( * (char**) supernode1, * (char**) supernode2);
}
