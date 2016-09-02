/*!
   \file simulation.c
   \brief Source code for simulation.h
   \author Piyush Divyankar
   \date 01/09/2016
 */

#include "simulation.h"

/** Boltzmann contant in eV/K */
#define KB  8.61e-5

void twoPhaseEquilibriaSimulation(unsigned long int seed_value)
{
    parameter * AlNi_fcc = _defaultFCCparameter();
    binEAMpot * potential = NULL;

    potential = eam_data_read("file_list.txt", "Al", "Ni");
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;

    const gsl_rng_type * T;
    gsl_rng * r;


    char * ouput_file_name = "AlNi_fcc-result.crystal.fcc";
    FILE * fp = fopen(ouput_file_name, "w");
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);

    for (size_t i = 0; i < steps; i++)
    {
        double u    = gsl_rng_uniform(r);
        int index   = u * AlNi_fcc->no_of_atoms;
        double e = energyToSwap(index, inputMatrix, potential,
                                AlNi_fcc, fccNeighbours);
        u = gsl_rng_uniform(r);
        e = exp(e/(KB * AlNi_fcc->temperature));
        if (e > 1 || e < u) {
            if (inputMatrix[index] == 0) {
                inputMatrix[index] = 1;
            }
            else inputMatrix[index] = 0;
        }
    }

    fwrite(inputMatrix, sizeof(inputMatrix), 1, fp);
}
