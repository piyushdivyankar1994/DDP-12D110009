/*!
   \file chemicalPotential.h
   \brief Function prototypes for chemicalPotential calculation
   \author Piyush Divyankar
   \date 1/08/2016
 */

#include "eam.h"

#ifndef _CHEM_POT_H
#define _CHEM_POT_H 1
double chemicalPotentialAtIndex(int, int *, binEAMpot *, parameter *, Sn_fcc *);
float concentrationOfAtom1(int *, int);
#endif
