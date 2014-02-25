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
  maths.c - cortex_var mathematical functions
*/

#include <stdlib.h>
#include <math.h>

// cortex_var headers
#include "maths.h"


// log(n!)= sum from i=1 to n, of  (log(i))                                                                                                                                                                          
float log_factorial(int number)
{
  if (number<0)
    {
      die("Do not call log_factorial with negative argument %d\n", number);
    }
  int i;
  float ret=0;
  for (i=1; i<=number; i++)
    {
      ret+=log(i);
    }
  return ret;
}


float log_factorial_ll(long long number)
{
  if (number<0)
    {
      die("Do not call log_factorial with negative argument %lld\n", number);
    }
  long long i;
  float ret=0;
  for (i=1; i<=number; i++)
    {
      ret+=log(i);
    }
  return ret;
}

int min_of_ints(int a, int b)
{
  if (a<b)    
    {
      return a;
    }
  else
    {
      return b;
    }
}

int max_of_ints(int a, int b)
{
  if (a>=b)    
    {
      return a;
    }
  else
    {
      return b;
    }
}

unsigned long calculate_mean_ulong(unsigned long* array, unsigned long len)
{
  unsigned long sum = 0;
  unsigned long num = 0;
  unsigned long i;

  for(i = 0; i < len; i++)
  {
     sum += i*array[i];
    num += array[i];
  }

  return num == 0 ? 0 : (sum / num);
}

long long calculate_mean(long long* array, long long len)
{
  long long sum=0;
  long long num=0;
  long long i;
  for (i=0; i<len; i++)
    {
      sum += i*array[i];
      num += array[i];
    }
  if (num>0)
    {
      return  (sum/num);
    }
  else
    {
      return 0;
    }
}

/*
void set_int_array_to_zero(int* array, int len)
{
  int i;
  for (i=0; i<len; i++)
    {
      array[i]=0;
    }
}
*/
