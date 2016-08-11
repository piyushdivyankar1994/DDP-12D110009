/*#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>*/
#include "eam_functions.h"
#include "math_functions.h"
#include "utils.h"
/*
 *follwing are the index fields and value feilds
*anything other than this is considered illegal and function
 */
//function declarations

//main function to wrap all the functions
//invoking main will follow by two things input file name and parameters file name
int neighbour_lattice_sites_number[]={12,6,24,12,24,8,48};

int main(int argc, char **argv)
{
    //making so that main takes input file name from bash this makes running test runs through bash script easier
    char *input_file,*parameter_file;
    argv++;
    input_file=*argv;
    argv++;
    parameter_file=*argv;


    double eam_data[45000];
   // int *site[4];

    double sites[ TOTAL_NEIGHBOUR_ATOMS_FCC_ARRAY_SIZE ];

    neighbour_lattice_sites_read(sites);

    eam_data_read(eam_data,"file_list.txt");       //to be replaced by new function and structure
    //printf("%d\n",neighbour_lattice_sites_number[3] );
    eam_monte_carlo_simulation(eam_data,sites,neighbour_lattice_sites_number,7,100000, parameter_file,input_file,"Result_lattice");

    return 0;
}
