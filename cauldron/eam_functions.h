#ifndef _EAM_FUNCTIONS_H
  #define  _EAM_FUNCTIONS_H
  #include<stdio.h>
  #include<stdlib.h>
  #include <gsl/gsl_rng.h>
  #include<math.h>
  #include "math_functions.h"
  #include "eam_functions.h"
  #include "utils.h"

  //obselete
  int field_check_eam(int ,int);
  //obselete use structure retrive functions instead which has better time complexity
  double eam_data_interpolation_func(double* , int , int , double );
  // broken up into other functions 
  void eam_monte_carlo_simulation(double *,double *,int *,int,int, char *, char *,char *);
  //
  void neighbour_lattice_sites_read(double* );
  // obselete replaced by a structure follow the structure is datatypes.h
  int read_parameter_file(char*,double []);

  void random_crystal_generator(char*, char* );
  // replaced with structure
  void eam_data_read(double*,char*);


  #define TOTAL_NEIGHBOUR_ATOMS_FCC 134
  #define TOTAL_NEIGHBOUR_ATOMS_FCC_ARRAY_SIZE 402
  #define NORMAL_TEMPERATURE 2000

#endif

#define RADIUS 0
#define PAIRING_POTENTIAL_Al_Al 1
#define PAIRING_POTENTIAL_Ni_Al 2
#define PAIRING_POTENTIAL_Ni_Ni 3
#define ELECTRON_DENSITY_Al 4
#define ELECTRON_DENSITY_Ni 5
#define ELECTRON_DENSITY_INDEX 6
#define EMBEDDING_FUNCTION_Ni 7
#define EMBEDDING_FUNCTION_Al 8

#define BOLTZMANN_CONST 1.38E-23
#define ELECTRONIC_CHARGE 1.602e-19

#define Al 1
#define Ni 2
