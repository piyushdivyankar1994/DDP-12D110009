/*!
   \file simulation.c
   \brief Source code for simulation.h
   \author Piyush Divyankar
   \date 01/09/2016
 */

#include "simulation.h"
#include "constants.h"
#include "gsl/gsl_statistics.h"
/** Boltzmann contant in eV/K */
#define KB 8.617e-5

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
    lookUpTable * t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
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
        int r = rand() % 12;
        point3D * current = point3D_indexToPoint3D_fcc(index, AlNi_fcc);
        point3D * ngbr = point3D_addVectors(current, &(fccNeighbours->s1n[r]));
        point3D_periodicBoundaryTransform(ngbr, AlNi_fcc);
        int ngbrIndex = point3D_point3DtoIndexFCC(ngbr, AlNi_fcc);
        inputMatrix[ngbrIndex] = 1 - inputMatrix[ngbrIndex];
        double e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);


        double p = exp(-(e2 - e1) / (KB * AlNi_fcc->temperature));
        if (p > 1 || p > u)
        {
            // Accept the swap
            swap_count++;
        }
        else
        {
            inputMatrix[index] = 1 - inputMatrix[index];
            inputMatrix[ngbrIndex] = 1 - inputMatrix[ngbrIndex];
        }
        free(current);
        free(ngbr);
        if (i % (AlNi_fcc->no_of_atoms / 100) == 0)
        {
            printf("%zu\t%le\n", i, totalEnergyQuick(inputMatrix, AlNi_fcc, fccNeighbours));
        }
    }
// /    printf("%f %f\n", concentration, (float)swap_count/(float)steps);
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

    lookUpTable * t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
    buildInstantEnergyLookup(t, potential);

    AlNi_fcc->N_MCS = 25;
    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
    print_parameters(AlNi_fcc);
    printf("#--------------------------------------------------------------\n");
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

    // int N1 = atomsType1(inputMatrix, AlNi_fcc);

    // Simulation constants

    /** Wave Number */
    double k = 0.2;
    double mu_0 = 0.2;

    double a = mu_0 + 2 * R * k * AlNi_fcc->temperature * (double) AlNi_fcc->no_of_atoms;
    double b = .05 * k * AlNi_fcc->temperature;

    double * chemPot = malloc(sizeof(double) * AlNi_fcc->no_of_atoms);
    for (int i = 0; i < AlNi_fcc->no_of_atoms; i++)
    {
        chemPot[i] = 0;
    }
    printf("## a = %0.3f\n ## b = %0.3f\n", a, b);
    int * ngbrTable = point3D_neighbourIndexTable_FCC(AlNi_fcc, fccNeighbours, THIRD);
    double * concTable = createConcentrationTable(inputMatrix, AlNi_fcc, ngbrTable, THIRD_NEAREST_FCC);
    for (size_t i = 0; i < steps; i++)
    {
        /* Step1: Selecting a random index */
        double u    = gsl_rng_uniform(r);
        int index   = u * AlNi_fcc->no_of_atoms;
        /* 1a. calculating energy  */
        //double e1 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);
        /* Step2: Flipping its spin */
        //inputMatrix[index] = 1 - inputMatrix[index];
        /* 2a. energy after switch */
        //double e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);

        double p;// = e2 - e1;
        /* We look at change in atoms of type 1 */

        if (inputMatrix[index] == 1)
        {
            p = (a + b * (84 * concTable[index] + 1));
        }
        else
        {
            p = (a + b * (84 * concTable[index] - 1));
        }

        // p = exp(-p / (KB * AlNi_fcc->temperature));
        printf("%le\n", p);
        u = gsl_rng_uniform(r);
        //if(u > 0.5)
        //inputMatrix[index] = 1 - inputMatrix[index];
        /*
        if (p > 1 || p > u)
        {
            if (inputMatrix[index] == 1)
            {
                updateConcentrationTable(concTable, ngbrTable, index, 42, 1);
            }
            else
            {
                updateConcentrationTable(concTable, ngbrTable, index, 42, 0);
            }
        }
        else
        {
            inputMatrix[index] = 1 - inputMatrix[index];
        }*/
    }
}

/**
 * Function lattice parameter simulations.
 */
void latticeParameterSimulation(size_t seed_value, double temp)
{
    parameter * AlNi_fcc = parameterReadFromFile("defaultFCC.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();
    AlNi_fcc->Nx = 4;
    AlNi_fcc->Ny = 4;
    AlNi_fcc->Nz = 4;
    AlNi_fcc->no_of_atoms = AlNi_fcc->Nx * AlNi_fcc->Ny * AlNi_fcc->Nz * 4;
    AlNi_fcc->N_MCS = 100;
    AlNi_fcc->temperature = temp;
    AlNi_fcc->lattice_parameter = 3.0;
    double * lc = malloc(sizeof(double) * AlNi_fcc->N_MCS);
    double * lc_std = malloc(sizeof(double) * AlNi_fcc->N_MCS);
    randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/lattice", 1729, 0);
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/lattice");
    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();
    // ----------------------------------------
    double stepSize = 0.00025;
    int count = 0;
    AlNi_fcc->temperature = 100;
    double pressure = 101325;
    double eOld = AlNi_fcc->no_of_atoms * energyAtIndexFCC(0, inputMatrix, potential, AlNi_fcc, fccNeighbours);
    double avgLP = 0;
    int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
    double * latticeP = malloc(sizeof(double) * AlNi_fcc->no_of_atoms);
    //print_parameters(AlNi_fcc);
    double eNew;
    int n = AlNi_fcc->Nx;

/*
    printf("#####################################################\n");
    printf("# MCS\t\ta\tvar_a\t\tstdDev_a\n" );
    printf("#####################################################\n");
*/
    for (int i = 0; i < steps; i++)
    {
        double u = gsl_rng_uniform(r);
        double da = (2 * u - 1) * stepSize;
        AlNi_fcc->lattice_parameter += da;
        double eNew = AlNi_fcc->no_of_atoms * energyAtIndexFCC(0, inputMatrix, potential, AlNi_fcc, fccNeighbours);
        double a = AlNi_fcc->lattice_parameter;
        double enthalpy = eNew - eOld + n*n*n*(pressure*3*a*a*da*6.242e-12 - 12 * KB * AlNi_fcc->temperature * log(1 + da/a));
        double p = exp(enthalpy / (-KB * AlNi_fcc->temperature)) > 1? 1 : exp(enthalpy / (-KB * AlNi_fcc->temperature));

        u = gsl_rng_uniform(r);
        if (u < p)
        {
            count++;
        }
        else{
            AlNi_fcc->lattice_parameter -= da;
        }
        //printf("%f\t%f\t%f\t%f\t%f\t\n", enthalpy, p, u, a, eNew-eOld);
        eOld = eNew;
        avgLP += AlNi_fcc->lattice_parameter;
        latticeP[i % AlNi_fcc->N_MCS] = AlNi_fcc->lattice_parameter;
        if (i % AlNi_fcc->no_of_atoms == 0 && i > 0)
        {
            lc[i / AlNi_fcc->no_of_atoms] = gsl_stats_mean(latticeP, 1, AlNi_fcc->N_MCS);
            lc_std[i / AlNi_fcc->no_of_atoms] = sqrt(gsl_stats_variance(latticeP, 1, AlNi_fcc->N_MCS));
            //printf("%d\t%f\t%f\t%f\n", i / AlNi_fcc->no_of_atoms, lc[i / AlNi_fcc->no_of_atoms], pow(lc_std[i / AlNi_fcc->no_of_atoms],2.0), lc_std[i / AlNi_fcc->no_of_atoms]);
        }
    }
    //printf("## Acceptance rate:%f\n", (float)count/(float)steps);
    double a[2] = {0};
    for(int i = 40; i <100; i++) {
        a[0] += lc[i];
    }
    a[0] = a[0]/(60.0);
    printf("%f\t%f\n", a[0] , temp);
}

void semiGrandCanonical_concentration_study(size_t seed_value)
{
    // Loading Resources
    parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

    AlNi_fcc->lattice_parameter = 3.65;
    lookUpTable * t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
    buildInstantEnergyLookup(t, potential);

    AlNi_fcc->N_MCS = 10;
    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
    print_parameters(AlNi_fcc);
    printf("#--------------------------------------------------------------\n");
    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();
    // ----------------------------------------

    double concentration = 0.04 ;   /// save plots for 0.5
    AlNi_fcc->temperature = 5;

    randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/input", rand(), concentration);
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");

    int N1 = atomsType1(inputMatrix, AlNi_fcc);     // number of of Al in the matrix
    float dc = 1.0/3.2e4;                           // change in concentration at each step this is absolute value
    printf("## %d, %le\n", N1, dc);

    int index;
    double e1, e2, ttl, p;
    double mu = -1e5;                 // this is mostly no of atoms times the chemical potential
    while(mu < 4.5e4) {
    double chemE = mu * dc;           // chemical driving force
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
    for (size_t i = 0; i < steps; i++) {
      double u = gsl_rng_uniform(r);
      index = u * AlNi_fcc->no_of_atoms;      // selecting random latice site
      e1 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours); // energy before flipping
      inputMatrix[index] = 1 - inputMatrix[index];          // flipping the spin at site
      e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);  // energy after flipping


      /**
       * If current state is 1 i.e after flipping atom became nickel we write
       * total hamiltonian as follows. Here we subtract chemical driving force.
       */
      if(inputMatrix[index] == 1) ttl = (e2 - e1) - chemE;


      else ttl = (e2 - e1) + chemE;

      ttl = exp(-(ttl)/(KB * AlNi_fcc->temperature));
      if( ttl < 1) {
        p = gsl_rng_uniform(r);
        if(p > ttl) {inputMatrix[index] = 1 - inputMatrix[index];}
      }
      //printf("%le\n", ttl);
    }
    N1 = atomsType1(inputMatrix, AlNi_fcc);
    printf("%f, %f\n", mu/32000.0, (float)N1/(float)AlNi_fcc->no_of_atoms);
    if(mu > 0) mu = mu + 1.25e3;
    else mu += 3e3;
  }

    double mu_0 = 0.2;
}

void pairwiseConstants() {
  parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
  binEAMpot * potential = NULL;

  potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
  Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

  AlNi_fcc->lattice_parameter = 3.65;
  double r = AlNi_fcc->lattice_parameter / sqrt(2.0);
  int k = binary_search(potential->radius, 0, 4999, r);
  printf("## Row position in Table:\t\t%d\n", k);
  printf("## E-Al-Al 1st nearest:\t\t%f\n", linear_interpolator(potential->radius[k],
                 potential->pair_atom1_atom1[k], potential->radius[k+1],
                 potential->pair_atom1_atom1[k+1], r));
  printf("## E-Ni-Al 1st nearest:\t\t%f\n", linear_interpolator(potential->radius[k], potential->pair_atom1_atom2[k], potential->radius[k+1], potential->pair_atom1_atom2[k+1], r));
  printf("## E-Ni-Ni 1st nearest:\t\t%f\n", linear_interpolator(potential->radius[k],
                 potential->pair_atom2_atom2[k], potential->radius[k+1],
                 potential->pair_atom2_atom2[k+1], r));
  r *= sqrt(2.0);
  k = binary_search(potential->radius, 0, 4999, r);
  printf("## Row position in Table:\t\t%d\n", k);
  printf("## E-Al-Al 2nd nearest:\t\t%f\n", linear_interpolator(potential->radius[k],
                 potential->pair_atom1_atom1[k], potential->radius[k+1],
                 potential->pair_atom1_atom1[k+1], r));
  printf("## E-Ni-Al 2nd nearest:\t\t%f\n", linear_interpolator(potential->radius[k], potential->pair_atom1_atom2[k], potential->radius[k+1], potential->pair_atom1_atom2[k+1], r));
  printf("## E-Ni-Ni 2nd nearest:\t\t%f\n", linear_interpolator(potential->radius[k],
                 potential->pair_atom2_atom2[k], potential->radius[k+1],
                 potential->pair_atom2_atom2[k+1], r));
}


void semiGrandCanonical_concentration_study_ljp(size_t seed_value)
{
    // Loading Resources
    parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

    AlNi_fcc->lattice_parameter = 3.65;
    lookUpTable * t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
    buildInstantEnergyLookup(t, potential);

    AlNi_fcc->N_MCS = 10;
    long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
    print_parameters(AlNi_fcc);
    printf("#--------------------------------------------------------------\n");
    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed_value);
    gsl_rng_env_setup();
    // ----------------------------------------

    double concentration = 0.1 ;   /// save plots for 0.5

    //NOTE: *****************************************
    AlNi_fcc->temperature = 5;//******************

    randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/input", rand(), concentration);
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
    // printf("## initial concentration = %f\n", (float)atomsType1(inputMatrix, AlNi_fcc)/(float)AlNi_fcc->no_of_atoms);
    const int N1 = atomsType1(inputMatrix, AlNi_fcc);     // number of of Al in the matrix
    float dc = 1.0/3.2e4;                           // change in concentration at each step this is absolute value
    // printf("## atoms of aluminum =%d\n", N1);

    int index;
    double e1, e2, ttl, p;
    double mu = -1e5;                 // this is mostly no of atoms times the chemical potential

    /** This is a quick access table the lists neighbours for all the sites in lattice */
    int * ngbrTable = point3D_neighbourIndexTable_FCC(AlNi_fcc, fccNeighbours, 2);

    /** Following are some constants specific to purely electrosatic type potentials */
    double E11_1 = -10;
    double E22_1 = -10;
    double E12_1 = -9.7;
    double E11_2 = -2;
    double E22_2 = -2;
    double E12_2 = -2;
    double ljp_temperature = 0.8617 / KB;
    int flag = 0;
    int Accepted = 0;

    while(concentration < 1) {
      Accepted = 0;
      randomMatrixGeneratorFCC(AlNi_fcc, "inputCrystalFiles/input", rand(), concentration);
      ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
      for (size_t i = 0; i < steps; i++) {
        double u;
        while(1) {
          u = gsl_rng_uniform(r);
          index = u * AlNi_fcc->no_of_atoms;      // selecting random latice site
          if(inputMatrix[index] == flag) {
            break;
          }
        }

        /*
        e1 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours); // energy before flipping
        inputMatrix[index] = 1 - inputMatrix[index];          // flipping the spin at site
        e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);  // energy after flipping
        */
        int s1=0, s2 =0;
        for(int j = 0; j < 18; j++){
          if(j < 12 && inputMatrix[ngbrTable[index*18+j]] == 0) s1++;
          else if(inputMatrix[ngbrTable[index*18+j]] == 0) s2++;
        }

        if(inputMatrix[index]==0){
          e1 = s1 * E11_1 + (12-s1) * E12_1 + s2 * E11_2 + (6-s2) * E12_2;
          e2 = s1 * E12_1 + (12-s1) * E22_1 + s2 * E12_2 + (6-s2) * E22_2;
        }
        else {
          e1 = s1 * E12_1 + (12-s1) * E22_1 + s2 * E12_2 + (6-s2) * E22_2;
          e2 = s1 * E11_1 + (12-s1) * E12_1 + s2 * E11_2 + (6-s2) * E12_2;
        }



        ttl = exp(-(e2-e1)/(KB * ljp_temperature));

        if( ttl < 1) {
          p = gsl_rng_uniform(r);
          if(ttl < p) {
            /** Here it means switch happened */
            inputMatrix[index] = 1 - inputMatrix[index];
            flag = 1 - flag;
            Accepted++;
          }
          /** if switch doesn't happen it means inputMatrix[index] remains same and flag remains same */
        }
        else {
          /** Here it means switch happened */
          inputMatrix[index] = 1 - inputMatrix[index];
          flag = 1 - flag;
          Accepted++;
        }
      }
      printf("%f %f\n", concentration, (float)Accepted/3200.0);
      concentration += 0.1;
    }
}

void ljp_cannonical_order_disorder_transformations(size_t seed_value) {
  // Loading Resources
  parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
  binEAMpot * potential = NULL;

  potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
  Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

  AlNi_fcc->lattice_parameter = 3.65;
  lookUpTable * t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
  buildInstantEnergyLookup(t, potential);

  AlNi_fcc->N_MCS = 10;
  long long int steps = AlNi_fcc->N_MCS * AlNi_fcc->no_of_atoms;
  print_parameters(AlNi_fcc);
  printf("#--------------------------------------------------------------\n");
  // Random number generator initialization
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed_value);
  gsl_rng_env_setup();
  // ----------------------------------------

  double concentration = 0.1 ;   /// save plots for 0.5

  //NOTE: *****************************************
  AlNi_fcc->temperature = 5;//******************

  L12lattice(AlNi_fcc, "inputCrystalFiles/input", rand());
  ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
  // printf("## initial concentration = %f\n", (float)atomsType1(inputMatrix, AlNi_fcc)/(float)AlNi_fcc->no_of_atoms);
  const int N1 = atomsType1(inputMatrix, AlNi_fcc);     // number of of Al in the matrix
  float dc = 1.0/3.2e4;                           // change in concentration at each step this is absolute value
  // printf("## atoms of aluminum =%d\n", N1);

  int index;
  double e1, e2, ttl, p;
  double mu = -1e5;                 // this is mostly no of atoms times the chemical potential

  /** This is a quick access table the lists neighbours for all the sites in lattice */
  int * ngbrTable = point3D_neighbourIndexTable_FCC(AlNi_fcc, fccNeighbours, 2);

  /** Following are some constants specific to purely electrosatic type potentials */
  double E11_1 = -10;
  double E22_1 = -10;
  double E12_1 = -9.7;
  double E11_2 = -2;
  double E22_2 = -2;
  double E12_2 = -2;
  double ljp_temperature = 0.1 / KB;
  int flag = 0;
  int Accepted = 0;

  while(ljp_temperature < 1 / KB) {
    Accepted = 0;
    L12lattice(AlNi_fcc, "inputCrystalFiles/input", rand());
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
    printf("%d\n", orderedPhaseCount(inputMatrix, AlNi_fcc, fccNeighbours));
    for (size_t i = 0; i < steps; i++) {
      double u;
      while(1) {
        u = gsl_rng_uniform(r);
        index = u * AlNi_fcc->no_of_atoms;      // selecting random latice site
        if(inputMatrix[index] == flag) {
          break;
        }
      }

      /*
      e1 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours); // energy before flipping
      inputMatrix[index] = 1 - inputMatrix[index];          // flipping the spin at site
      e2 = energyAtIndexFCC_fast(index, inputMatrix, AlNi_fcc, fccNeighbours);  // energy after flipping
      */
      int s1=0, s2 =0;
      for(int j = 0; j < 18; j++){
        if(j < 12 && inputMatrix[ngbrTable[index*18+j]] == 0) s1++;
        else if(inputMatrix[ngbrTable[index*18+j]] == 0) s2++;
      }

      if(inputMatrix[index]==0){
        e1 = s1 * E11_1 + (12-s1) * E12_1 + s2 * E11_2 + (6-s2) * E12_2;
        e2 = s1 * E12_1 + (12-s1) * E22_1 + s2 * E12_2 + (6-s2) * E22_2;
      }
      else {
        e1 = s1 * E12_1 + (12-s1) * E22_1 + s2 * E12_2 + (6-s2) * E22_2;
        e2 = s1 * E11_1 + (12-s1) * E12_1 + s2 * E11_2 + (6-s2) * E12_2;
      }



      ttl = exp(-(e2-e1)/(KB * ljp_temperature));

      if( ttl < 1) {
        p = gsl_rng_uniform(r);
        if(ttl < p) {
          /** Here it means switch happened */
          inputMatrix[index] = 1 - inputMatrix[index];
          flag = 1 - flag;
          Accepted++;
        }
        /** if switch doesn't happen it means inputMatrix[index] remains same and flag remains same */
      }
      else {
        /** Here it means switch happened */
        inputMatrix[index] = 1 - inputMatrix[index];
        flag = 1 - flag;
        Accepted++;
      }
    }
    printf("%f, %d \n", ljp_temperature * KB ,orderedPhaseCount(inputMatrix, AlNi_fcc, fccNeighbours));
    ljp_temperature += 0.1/KB;
  }
}

void bccCannonicalBenchmark() {
    // Loading Resources
    parameter * pBcc = _defaultBCCparameter();

    Sn_bcc * bccNgbrs = readBCCfromFile( "/home/piyush/Desktop/DDP-12D110009/neighbours/bcc/bccNeighbours.txt");

    pBcc->N_MCS = 2000;
    long long int steps = pBcc->N_MCS * pBcc->no_of_atoms;
    print_parameters(pBcc);
    printf("#--------------------------------------------------------------\n");
    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 465653);
    gsl_rng_env_setup();
    // ----------------------------------------

    randomMatrixGeneratorFCC(pBcc, "inputCrystalFiles/input", rand(), 0.5   );
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");

    double E11_1 = -10;
    double E22_1 = -10;
    double E12_1 = -9.7;
    double E11_2 = -2;
    double E22_2 = -2;
    double E12_2 = -2;
    double de1, de2;
    printf("## E11_1 = %f\t\n", E11_1);
    printf("## E12_1 = %f\t\n", E12_1);
    printf("## E22_1 = %f\t\n", E22_1);
    printf("## E11_2 = %f\t\n", E11_2);
    printf("## E12_2 = %f\t\n", E12_2);
    printf("## E22_2 = %f\t\n", E22_2);
    double ljp_temp = 0.817;
    int * ngbrTable = point3D_neighbourIndexTable_BCC(pBcc, bccNgbrs);
    int Accepted1 = 0;
    int Accepted2 = 0;
    double u;
    int index1, index2;
    //int flag = 0;
    double  a1 = 0, b1 = 0, a2 = 0, b2 = 0;
    printf("## concentration = %f\n", concentrationOfAtom1(inputMatrix, pBcc->no_of_atoms));
    for (size_t i = 0; i < steps; i++)
    {
        /* Step1: Selecting a random index */
        index1 = gsl_rng_uniform(r) * pBcc->no_of_atoms;
        while(1) {
            u = gsl_rng_uniform(r);
            index2 = u * pBcc->no_of_atoms;      // selecting random latice site
            if(inputMatrix[index1] != inputMatrix[index2]) break;   //with opposite spin
          }

        a1 = 0;
        b1 = 0;
        a2 = 0;
        b2 = 0;
        for(int j = 0; j < 8; j++) {
            if(inputMatrix[ngbrTable[index1*14 + j]] == 0) a1++;
            if(inputMatrix[ngbrTable[index2*14 + j]] == 0) b1++;
        }
        for(int j = 8; j < 14; j++) {
            if(inputMatrix[ngbrTable[index1*14 + j]] == 0) a2++;
            if(inputMatrix[ngbrTable[index1*14 + j]] == 0) b2++;
        }
        de1 = (a1*(E12_1 - E11_1) + (8-a1)*(E22_1 - E12_1) + a2*(E12_2 - E11_2) \
                + (6-a2)*(E22_2 - E12_2))*0.5;
        de2 = (b1*(E11_1 - E12_1) + (8-b1)*(E12_1 - E22_1) + b2*(E11_2 - E12_2) \
                + (6-b2)*(E12_2 - E22_2))*0.5;

        if(inputMatrix[index1] == 1) {
            de1 = -de1;
            de2 = -de2;
        }
        double ar = exp(-(de1+de2)/(ljp_temp));

        if( ar < 1) {
          double p = gsl_rng_uniform(r);
          if(ar < p) {
            /** Here it means switch happened */
            inputMatrix[index1] = 1 - inputMatrix[index1];
            inputMatrix[index2] = 1 - inputMatrix[index2];
            if(i/pBcc->no_of_atoms > 1000)Accepted1++;
          }
          /** if switch doesn't happen it means inputMatrix[index] remains same and flag remains same */
        }
        else {
          /** Here it means switch happened */
          inputMatrix[index1] = 1 - inputMatrix[index1];
          inputMatrix[index2] = 1 - inputMatrix[index2];
          if(i/pBcc->no_of_atoms > 1000)Accepted2++;
      }
    }
    printf("## Acceptance rate at 0.5 ==%f,%f\n", ((float)Accepted1/(float)steps)*200, ((float)Accepted2/(float)steps)*200 );
    printf("## concentration = %f\n", concentrationOfAtom1(inputMatrix, pBcc->no_of_atoms));
}

void bcc_SGCannonicalBenchmark() {
    // Loading Resources
    parameter * pBcc = _defaultBCCparameter();

    Sn_bcc * bccNgbrs = readBCCfromFile( "/home/piyush/Desktop/DDP-12D110009/neighbours/bcc/bccNeighbours.txt");

    pBcc->N_MCS = 20000;
    long long int steps = pBcc->N_MCS * pBcc->no_of_atoms;
    print_parameters(pBcc);
    printf("#--------------------------------------------------------------\n");
    // Random number generator initialization
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 10000);
    gsl_rng_env_setup();
    // ----------------------------------------

    double E11_1 = -10;
    double E22_1 = -10;
    double E12_1 = -9.7;
    double E11_2 = -2;
    double E22_2 = -2;
    double E12_2 = -2;
    double de1, de2;
    printf("## E11_1 = %f\t\n", E11_1);
    printf("## E12_1 = %f\t\n", E12_1);
    printf("## E22_1 = %f\t\n", E22_1);
    printf("## E11_2 = %f\t\n", E11_2);
    printf("## E12_2 = %f\t\n", E12_2);
    printf("## E22_2 = %f\t\n", E22_2);
    double ljp_temp = 0.817;
    int * ngbrTable = point3D_neighbourIndexTable_BCC(pBcc, bccNgbrs);
    int Accepted1 = 0;
    int Accepted2 = 0;
    double u;
    int index1, index2;
    //int flag = 0;
    double  a1 = 0, b1 = 0, a2 = 0, b2 = 0;
    double dmu = -0.8;
    printf("## dmu, Accepted, cAl");
    while(dmu < 0.9) {
        Accepted2 = 0;
        Accepted1 = 0;
        randomMatrixGeneratorFCC(pBcc, "inputCrystalFiles/input", rand(), 0.5   );
        ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
        for (size_t i = 0; i < steps; i++)
        {
            /* Step1: Selecting a random index */
            index1 = gsl_rng_uniform(r) * pBcc->no_of_atoms;

            a1 = 0;
            a2 = 0;
            for(int j = 0; j < 8; j++) {
                if(inputMatrix[ngbrTable[index1*14 + j]] == 0) a1++;
            }
            for(int j = 8; j < 14; j++) {
                if(inputMatrix[ngbrTable[index1*14 + j]] == 0) a2++;
            }
            de1 = (a1*(E12_1 - E11_1) + (8-a1)*(E22_1 - E12_1) + a2*(E12_2 - E11_2) \
                    + (6-a2)*(E22_2 - E12_2))*0.5 + dmu;

            if(inputMatrix[index1] == 1) {
                de1 = -de1;
            }
            double ar = exp(-(de1)/(ljp_temp));

            if( ar < 1) {
              double p = gsl_rng_uniform(r);
              if(ar < p) {
                /** Here it means switch happened */
                inputMatrix[index1] = 1 - inputMatrix[index1];
                Accepted1++;
              }
              /** if switch doesn't happen it means inputMatrix[index] remains same and flag remains same */
            }
            else {
              /** Here it means switch happened */
              inputMatrix[index1] = 1 - inputMatrix[index1];
              Accepted2++;
          }
        }
        printf("%f, %f, %f\n", dmu, ((float)(Accepted1+Accepted2)/(float)steps)*100, concentrationOfAtom1(inputMatrix, pBcc->no_of_atoms));
        dmu += 0.1;
    }
}
