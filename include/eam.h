/*!
   \file eam.h
   \brief Function prototypes of EAM methods
   \author Piyush Divyankar
   \date 1/08/2016
 */

#include "datatypes.h"
#include "point3D.h"

#ifndef _EAM_H
  #define _EAM_H 1

int * atomicMatrixRead(char *, parameter *);
int * readCrystalFileFCC(char *);

void print_AtomicMatrix(int *, int, int);

Sn_fcc * Sn_fcc_readNeighbours_fromFile(char *);

void print_Neighbours(Sn_fcc *);

Sn_fcc * _defaultFCCNeighbours();

binEAMpot* eam_data_read(char *, char [], char []);

rdf * rdf_radius_retrive(binEAMpot *, double);

eDen_df * eDen_df_charge_density_retrive(binEAMpot *, double);

double energyAtIndexFCC(int, int *, binEAMpot *, parameter *, Sn_fcc *);

double* energyInMatrix(int *, binEAMpot *, parameter *, Sn_fcc *);

void printEnergyMap(double *, int, int);

double energyToSwap(int, int *, binEAMpot *, parameter *, Sn_fcc *);
double energyToSwap_fast(int index, int * a, parameter * p, Sn_fcc * ngbrs);

double* deltaEnergyMatrix(int *, binEAMpot *, parameter *, Sn_fcc *);

lookUpTable* createLookUpTable(binEAMpot *data, parameter *p, Sn_fcc *ngbrs);
void printLookUpTable(lookUpTable *p);

/** Temporary function no related to any code */
double * energyInMatrix_ver2(int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs, int noNgbrs);
/** Temporary function no related to any code */
double energyAtIndexFCC_ver2(int index, int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs, int noNgbrs);

double energyAtIndexFCC_fast(int, int *, parameter *, Sn_fcc *);
void buildInstantEnergyLookup(lookUpTable *data, binEAMpot *pot);

double avgConcentrationAtom1(ATOM *test, parameter *p);
int atomsType1(ATOM *test, parameter *p);

Sn_bcc * readBCCfromFile(char * fileName);


#endif /* ifndef _EAM_H */
