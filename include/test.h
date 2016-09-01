#include "datatypes.h"
#include "parameters.h"
#include "point3D.h"
#include "eam.h"
#include "chemicalPotential.h"
#include "analysis.h"

#ifndef _TEST_H
#define _TEST_H 1
  void test_new_parameters();
  void test_dataRetrival();
  void test_eam_data_read();
  void test_point3D_indexToPoint3D_fcc();
  void test_point3D();
  void test_defaultFCCparameter();
  void test_Sn_fcc_readNeighbours_fromFile();
  void test_AtomicMatrixRead();
  void test_energyAtIndexFCC();
  void test_point3D_point3DtoIndex();
  void test_point3D_periodicBoundaryTransform();
  void test_energyInMatrix();
  void test_deltaEnergyMatrix();
  void test_chemicalPotentialAtIndex();
  void test_analysis_totalEnergy();
  void test_orderedPhaseCount();
  void test_antiOrderedPhaseCount();
  void test_getRandom();
  void test_randomMatrixGeneratorFCC();
  void test_readCrystalFileFCC();

  void test_point3D_origin();
  void test_point3D_newPoint();
  void test_point3D_addVectors();

  void test_parametersInputOutput();
#endif
