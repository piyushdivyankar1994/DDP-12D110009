/*!
   \file simulation.c
   \brief Source code for simulation.h
   \author Piyush Divyankar
   \date 01/09/2016
 */

#include "simulation.h"
#include "gsl/gsl_statistics.h"
/** Boltzmann contant in eV/K */
#define KB  8.61e-5

/** Reservoir parameter */
const double R = 0;
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

    AlNi_fcc->N_MCS = 10;
    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
    AlNi_fcc->lattice_parameter = 3.56;
    double concentration = .2;

    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();

    randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/input", rand(), concentration);
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");

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
            // Accept the swap
            swap_count++;
        }
        else {
            inputMatrix[index] = 1 - inputMatrix[index];
            inputMatrix[ngbrIndex] = 1 - inputMatrix[ngbrIndex];
        }
        free(current);
        free(ngbr);
        if (i % (AlNi_fcc->no_of_atoms / 100) == 0) {
            printf("%zu\t%le\n", i, totalEnergyQuick(inputMatrix, AlNi_fcc, fccNeighbours));
        }
    }
///    printf("%f %f\n", concentration, (float)swap_count/(float)steps);
    free(inputMatrix);

    free(fccNeighbours);
    free(potential);
    free(AlNi_fcc);

    gsl_rng_free(r);
}


void semiGrandCanonical(size_t seed_value)
{
    // Loading Resources
    parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

    lookUpTable *t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
    buildInstantEnergyLookup(t, potential);

    AlNi_fcc->N_MCS = 1000;
    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;

    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();
    // ----------------------------------------

    double concentration = 0.5;

    randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/input", rand(), concentration);
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");

    int N1 = atomsType1(inputMatrix, AlNi_fcc);

    // Simulation constants

    /** Wave Number */
    double k = 0.2;
    double mu_0 = 0.2;

    double a = mu_0 + 2 * R * k * AlNi_fcc->temperature * (double)AlNi_fcc->no_of_atoms;
    double b = R * k * AlNi_fcc->temperature;

    for (size_t i = 0; i < steps; i++)
    {
        /* Step1: Selecting a random index */
        double u    = gsl_rng_uniform(r);
        int index   = u * AlNi_fcc->no_of_atoms;
        /* 1a. calculating energy  */
        double e1 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);
        /* Step2: Flipping its spin */
        inputMatrix[index] = 1 - inputMatrix[index];
        /* 2a. energy after switch */
        double e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);

        double p = e2 - e1;
        /* We look at change in atoms of type 1 */
        if (inputMatrix[index] == 1) {
            p = p - (a - b * (2 * N1 - 1))*(-0.5);
            N1--;
        }
        else {
            p = p - (a - b * (2 * N1 + 1))*(0.5);
            N1++;
        }

        p = exp(-p / (KB * AlNi_fcc->temperature));
        u = gsl_rng_uniform(r);
        if (p > 1 || p > u)
        {
            /* Accept the swap */
        }
        else
            inputMatrix[index] = 1 - inputMatrix[index];
    }
}

/**
 * Function generates and outputs lattice parameter v/s temperature plot.
 */
void latticeParameterSimulation(size_t seed_value) {
    parameter * AlNi_fcc = parameterReadFromFile("defaultFCC.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();
    AlNi_fcc->Nx = 4;
    AlNi_fcc->Ny = 4;
    AlNi_fcc->Nz = 4;
    AlNi_fcc->no_of_atoms = 64 * 4;
    AlNi_fcc->N_MCS = 10;
    AlNi_fcc->temperature = 100;

    randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/lattice", rand(), 0);
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/lattice");
    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();
    // ----------------------------------------
    double stepSize = 0.03;
    int count = 0;
    AlNi_fcc->temperature = 100;
    double pressure = 1e5;
    double eOld = AlNi_fcc->no_of_atoms * energyAtIndexFCC(0, inputMatrix, potential, AlNi_fcc, fccNeighbours);
    double avgLP = 0;
    int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
    double latticeP[AlNi_fcc->no_of_atoms];
    print_parameters(AlNi_fcc);

    for (int i = 0; i < steps; i++)
    {
        double u = gsl_rng_uniform(r);
        double da = (2 * u - 1) * stepSize;
        AlNi_fcc->lattice_parameter -= da;
        double eNew = AlNi_fcc->no_of_atoms * energyAtIndexFCC(0, inputMatrix, potential, AlNi_fcc, fccNeighbours);
        double enthalpy = eNew - eOld - KB * AlNi_fcc->temperature * 3 * AlNi_fcc->no_of_atoms *          \
                            log(1 + 3 * fabs(da) / AlNi_fcc->lattice_parameter) +    \
                            pressure * 3 * AlNi_fcc->lattice_parameter * AlNi_fcc->lattice_parameter * -da * 6.241e-12;

        double p = exp(enthalpy / (-KB * AlNi_fcc->temperature));

        if ((u < p && p < 1) || p > 1) {
            //Accepted
            count++;
        }
        else {
            AlNi_fcc->lattice_parameter += da;
        }
        eOld = eNew;
        avgLP += AlNi_fcc->lattice_parameter;
        latticeP[i % AlNi_fcc->N_MCS] = AlNi_fcc->lattice_parameter;
        if (i % AlNi_fcc->N_MCS == 0 && i > 0) {
            printf("%f,%d,%f\n",gsl_stats_mean(latticeP, 1, AlNi_fcc->N_MCS), i/AlNi_fcc->N_MCS,
                              gsl_stats_variance(latticeP, 1, AlNi_fcc->N_MCS));
        }
    }
    /*
    pressure = 1e5;
    double eOld = AlNi_fcc->no_of_atoms * energyAtIndexFCC(0, inputMatrix, potential, AlNi_fcc, fccNeighbours);
    double avgLP = 0;
    for(int i= 0; i < 1000; i += 10) {
        int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
        AlNi_fcc->lattice_parameter = 2.5;
        AlNi_fcc->N_MCS += 50;
        for (int i = 0; i < steps; i++)
        {
            double u = gsl_rng_uniform(r);
            double da = (2 * u - 1) * stepSize;
            AlNi_fcc->lattice_parameter -= da;
            double eNew = AlNi_fcc->no_of_atoms * energyAtIndexFCC(0, inputMatrix, potential, AlNi_fcc, fccNeighbours);
            double enthalpy = eNew - eOld - KB * AlNi_fcc->temperature * 3 * AlNi_fcc->no_of_atoms *          \
                                log(1 + 3 * fabs(da) / AlNi_fcc->lattice_parameter) +    \
                                pressure * 3 * AlNi_fcc->lattice_parameter * AlNi_fcc->lattice_parameter * -da * 6.241e-12;

            double p = exp(enthalpy / (-KB * AlNi_fcc->temperature));

            if ((u < p && p < 1) || p > 1) {
                //Accepted
                count++;
            }
            else {
                AlNi_fcc->lattice_parameter += da;
            }
            eOld = eNew;
            avgLP += AlNi_fcc->lattice_parameter;
            if (i % AlNi_fcc->no_of_atoms == 0 ) {
                latticeP[i / AlNi_fcc->no_of_atoms] = AlNi_fcc->lattice_parameter;
            }
        }
        printf("%d\t%f\t%f\n", AlNi_fcc->N_MCS, mean, sqrt(variance));
    }*/


}
