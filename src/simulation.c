#include "simulation.h"

void twoPhaseEquilibriaSimulation(unsigned long int seed_value) {
  parameter *AlNi_fcc = _defaultFCCparameter();
  binEAMpot *potential = NULL;
  eam_data_read(&(potential), "file_list.txt", "Al", "Ni");
  randomMatrixGeneratorFCC(AlNi_fcc, "test", seed_value + 213);
  
}
