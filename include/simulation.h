#include "datatypes.h"
#include "eam.h"
#include "point3D.h"
#include "analysis.h"
#include "parameters.h"
#include "math_functions.h"
#include "chemicalPotential.h"

#ifndef _SIMULATION_H
#define _SIMULATION_H 1
void cannonicalEnsemble(unsigned long int);
void semiGrandCanonical(size_t seed_value);
void latticeParameterSimulation(size_t seed_value, double temp);
void semiGrandCanonical_concentration_study(size_t seed_value);
void pairwiseConstants();

void ljp_cannonical_order_disorder_transformations();
void bccCannonicalBenchmark();
void bcc_SGCannonicalBenchmark();
#endif
