/*!
   \file simulation.c
   \brief Source code for simulation.h
   \author Piyush Divyankar
   \date 01/09/2016
 */

#include "simulation.h"

/** Boltzmann contant in eV/K */
#define KB  8.61e-5

/**
 * Simulates a cannoical ensemble.
 * @param seed_value Random seed value for simulation
 * @callgraph
 */
void cannonicalEnsemble(unsigned long int seed_value)
{
    parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

    lookUpTable *t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
    buildInstantEnergyLookup(t, potential);

    AlNi_fcc->N_MCS = 1000;
    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;

    double concentration = 0.1;

    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();


    while (concentration < 1.0)
    {
        randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/input", rand(), concentration);
        ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
        concentration += 0.1;

        int swap_count = 0;
        for (size_t i = 0; i < steps; i++)
        {
            double u    = gsl_rng_uniform(r);
            int index   = u * AlNi_fcc->no_of_atoms;
            u = gsl_rng_uniform(r);
            double e1 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);
            inputMatrix[index] = 1 - inputMatrix[index];
            int r = rand()%12;
            point3D *current = point3D_indexToPoint3D_fcc(index, AlNi_fcc);
            point3D *ngbr = point3D_addVectors(current, &(fccNeighbours->s1n[r]));
            point3D_periodicBoundaryTransform(ngbr, AlNi_fcc);
            int ngbrIndex = point3D_point3DtoIndexFCC(ngbr, AlNi_fcc);
            inputMatrix[ngbrIndex] = 1 - inputMatrix[ngbrIndex];
            double e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);


            double p = exp(-(e2 - e1)/(KB * AlNi_fcc->temperature));
            if (p > 1 || p > u)
            {
                /* Accept the swap */
                swap_count++;
            }
            else {
                inputMatrix[index] = 1 - inputMatrix[index];
                inputMatrix[ngbrIndex] = 1 - inputMatrix[ngbrIndex];
            }
            free(current);
            free(ngbr);
        }
        printf("%f %f\n", concentration, (float)swap_count/(float)steps);
        free(inputMatrix);
    }

    free(fccNeighbours);
    free(potential);
    free(AlNi_fcc);

    gsl_rng_free(r);
}
