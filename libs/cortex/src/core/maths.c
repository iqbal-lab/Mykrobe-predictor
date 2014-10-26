/*
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@tgac.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * The MIT License (MIT)
 * Copyright (c) 2009-2014 <Z. Iqbal and M. Caccamo>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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

uint64_t calculate_mean_uint64_t(uint64_t* array, uint64_t len)
{
  uint64_t sum = 0;
  uint64_t num = 0;
  uint64_t i;

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
