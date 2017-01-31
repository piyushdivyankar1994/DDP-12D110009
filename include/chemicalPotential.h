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
double concentrationAtIndex(ATOM * atomicMatrix, int * ngbrTable, int index, int n);
double * createConcentrationTable(ATOM * atomicMatrix, parameter * p, int * ngbrTable, int n);
void printConcentrationTable(double * table, int initial_index, int final_index);
void updateConcentrationTable(double * table, int * ngbrIndex, int index, int n, int change);
#endif
