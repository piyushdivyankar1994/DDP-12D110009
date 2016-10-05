/*!
   \file analysis.h
   \brief Function prototypes of analysis
   \author Piyush Divyankar
   \date 1/08/2016
 */
#include "eam.h"

#ifndef _ANALYSIS_H
#define _ANALYSIS_H 1
double analysis_totalEnergy(int *, binEAMpot *, parameter *, Sn_fcc *);
int orderedPhaseCount(int *, parameter *, Sn_fcc *);
int antiOrderedPhaseCount(int *, parameter *, Sn_fcc *);
void randomMatrixGeneratorFCC(parameter *, char *, unsigned long int, double);
int totalAtomsInFile(FILE *);
double totalEnergyQuick(ATOM * a, parameter * p, Sn_fcc * ngbrs);
#endif
