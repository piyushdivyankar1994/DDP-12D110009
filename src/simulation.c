/*!
   \file simulation.c
   \brief Source code for simulation.h
   \author Piyush Divyankar
   \date 01/09/2016
*/

#include "simulation.h"

void twoPhaseEquilibriaSimulation(unsigned long int seed_value)
{
    parameter * AlNi_fcc = _defaultFCCparameter();
    binEAMpot * potential = NULL;

    potential = eam_data_read("file_list.txt", "Al", "Ni");
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");

}
