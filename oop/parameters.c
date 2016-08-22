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
/// TODO: Make this more elegant
void createParameterFileFromInput() {
  parameter* new = (parameter*)malloc(sizeof(parameter));
  new = _defaultFCCparameter();
  printf("\nEnter # of MonteCarlo simulations per atom = ");
  scanf("%d", &(new->N_MCS));
  printf("\nEnter Lattice size in x-direction = ");
  scanf("%d", &(new->Nx));
  printf("\nEnter Lattice size in y-direction = ");
  scanf("%d", &(new->Ny));
  printf("\nEnter Lattice size in z-direction = ");
  scanf("%d", &(new->Nz));
  printf("\nEnter atoms per site = ");
  scanf("%d", &(new->atoms_per_site));
  printf("\nEnter lattice parameter = ");
  scanf("%le", &(new->lattice_parameter));
  printf("\nEnter simulation Temperature");
  scanf("%le", &(new->temperature));
  new->no_of_atoms = new->Nx * new->Ny * new->Nz * new->no_of_atoms;

  new->nearestNeighbours[0] = 12;
	new->nearestNeighbours[1] = 6;
	new->nearestNeighbours[2] = 24;
	new->nearestNeighbours[3] = 12;
	new->nearestNeighbours[4] = 24;
	new->nearestNeighbours[5] = 8;
	new->nearestNeighbours[6] = 48;
  char str[1000];
  printf("\nEnter a filename for list of parameters");
  scanf("%s\n", str);
  strcat(str, ".parameters");
  FILE *fp = fopen(str, "w");
  fwrite(new, sizeof(parameter), 1, fp);
  fclose(fp);
}
