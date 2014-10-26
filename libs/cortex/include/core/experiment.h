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
  experiment.h
*/

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include "global.h"

typedef enum
  {
    Unspecified                                                = 0,
    EachColourADiploidSample                                   = 1,
    EachColourADiploidSampleExceptTheRefColour                 = 2,
    EachColourAHaploidSample                                   = 3,
    EachColourAHaploidSampleExceptTheRefColour                 = 4,

    /*
    EachColourAPolyploidSample                                 = 5,
    EachColourPooledHaploidSamplesWithKnownNumberOfSamples     = 7,
    EachColourPooledDiploidSamplesWithKnownNumberOfSamples     = 7,
    EachColourPooledHaploidSamplesWithUnknownNumberOfSamples   = 9,
    EachColourPooledDiploidSamplesWithUnknownNumberOfSamples   = 9,
    MetaGenomic                                                = 10,*/
  } ExperimentType;

typedef enum {
  AssumeUncleaned                                = 1,
  AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice = 2,//because we cleaned off stuff that occurred only once
  NoIdeaWhatCleaning                             =3
}AssumptionsOnGraphCleaning;

#endif
