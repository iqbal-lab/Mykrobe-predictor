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
