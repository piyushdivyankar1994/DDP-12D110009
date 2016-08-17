#include "parameters.h"

void print_parameters(parameter* a){
  printf("Size in x-direction           = %d\n", a->Nx);
  printf("Size in y-direction           = %d\n", a->Ny);
  printf("Size in z-direction           = %d\n", a->Nz);
  printf("Lattice Parameter             = %le\n", a->lattice_parameter);
  printf("No. of MonteCarlo Simulations = %d\n", a->N_MCS);
  printf("Temperature                   = %le\n", a->temperature);
  printf("Atoms per site                = %d\n", a->atoms_per_site);
  printf("Totol No of atoms             = %d\n", a->no_of_atoms );
}


parameter* new_parameters(char* filename) {
  parameter* new = (parameter*)malloc(sizeof(parameter));
  FILE *fp;
  fp = fopen(filename, "r");
  int N;

  fscanf(fp, "%d", &N);
  new->Nx = N;
  printf("%d\n", N);

  fscanf(fp, "%d", &(new->Ny));
  fscanf(fp, "%d", &(new->Nz));
  fscanf(fp, "%le", &(new->lattice_parameter));
  fscanf(fp, "%d", &(new->N_MCS));
  fscanf(fp, "%le", &(new->temperature));
  fscanf(fp, "%d", &(new->atoms_per_site));

  new->nearestNeighbours[0] = 12;
  new->nearestNeighbours[1] = 6;
  new->nearestNeighbours[2] = 24;
  new->nearestNeighbours[3] = 12;
  new->nearestNeighbours[4] = 24;
  new->nearestNeighbours[5] = 8;
  new->nearestNeighbours[6] = 48;
  new->no_of_atoms = new->Nx * new->Ny * new->Nz * new->atoms_per_site;
  return new;
}

parameter* _defaultFCCparameter(){
	parameter* new = (parameter*)malloc(sizeof(parameter));
	new->N_MCS = 100;
	new->Nx	   = 10;
	new->Ny	   = 10;
	new->Nz	   = 10;
	new->atoms_per_site	   = 4;
	new->lattice_parameter = 4.00;
	new->temperature	   = 1000;
	new->nearestNeighbours[0] = 12;
	new->nearestNeighbours[1] = 6;
	new->nearestNeighbours[2] = 24;
	new->nearestNeighbours[3] = 12;
	new->nearestNeighbours[4] = 24;
	new->nearestNeighbours[5] = 8;
	new->nearestNeighbours[6] = 48;
	new->no_of_atoms 		  = new->Nx * new->Ny * new->Nz * new->atoms_per_site;
	return new;
}
