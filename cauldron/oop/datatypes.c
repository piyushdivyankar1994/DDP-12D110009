#include<stdio.h>
#include<stdlib.h>

struct parameter {
  int no_of_atoms;
  int Nx;
  int Ny;
  int Nz;
  double lattice_parameter;
  int N_MCS;
  double temperature;
  int atoms_per_site;
}typedef parameter;

struct lattice {
  float *atomic_positions;
  float *neighbour_coordinates;
  int *site_array;
}typedef lattice;

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
  double t;
  fscanf(fp, "%d", &N);
  new->Nx = N;
  printf("%d\n", N);

  fscanf(fp, "%d", &(new->Ny));
  fscanf(fp, "%d", &(new->Nz));
  fscanf(fp, "%le", &(new->lattice_parameter));
  fscanf(fp, "%d", &(new->N_MCS));
  fscanf(fp, "%le", &(new->temperature));
  fscanf(fp, "%d", &(new->atoms_per_site));


  new->no_of_atoms = new->Nx * new->Ny * new->Nz * new->atoms_per_site;
  return new;
}

void new_parameters_test(){
  parameter* new = NULL;
  new = new_parameters("parameter1.txt");
  print_parameters(new);
  free(new);
}

int main() {
  new_parameters_test();
}
