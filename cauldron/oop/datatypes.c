#include<stdio.h>
#include<stdlib.h>
#include "math_functions.h"

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


struct binaryEAMpotential {
  char atom1[2];
  int no_of_files;
  char atom2[2];
  double radius[5000];
  double atom1_charge_Density[5000];
  double atom2_charge_Density[5000];
  double pair_atom1_atom1[5000];
  double pair_atom1_atom2[5000];
  double pair_atom2_atom2[5000];
  double atom1_embedding_energy[5000];
  double atom2_embedding_energy[5000];
  double chargeDensity[5000];
  double minRadius;
  double maxRadius;
  double min_eDen;
  double max_eDen;
}typedef binEAMpot;

void eam_data_read(binEAMpot** eam_data, char *fileName, char atom1[2], char atom2[2])
{
    FILE *fp_init,*fp1;
    fp_init=fopen(fileName,"r");
    *eam_data = (binEAMpot*)malloc(sizeof(binEAMpot));
    //double* eam_data;
    //(*eam_dat)->atom1 = (char*)malloc(sizeof(atom1)+1);
    //(*eam_dat)->atom2 = (char*)malloc(sizeof(atom2)+1);
    (*eam_data)->atom1[0] = atom1[0];
    (*eam_data)->atom1[1] = atom1[1];

    (*eam_data)->atom2[0] = atom2[0];
    (*eam_data)->atom2[1] = atom2[1];
    (*eam_data)->no_of_files = 7;
    char str[100];
    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    double val1,val2;
    int j;
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j=0;j<5000;j++) {
      fscanf(fp1,"%le %le ",&val1,&val2);
      (*eam_data)->radius[j]        = val1;
      (*eam_data)->pair_atom1_atom1[j] = val2;
    }
    (*eam_data)->minRadius = (*eam_data)->radius[0];
    (*eam_data)->maxRadius = (*eam_data)->radius[4999];

    //printf("%d %le\n", j, (*eam_data)->radius[234]);
    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j = 0; j < 5000; j++) {
      fscanf(fp1, "%le %le ", &val1, &val2);
      //(*eam_data)->radius[j]        = val1;
      (*eam_data)->pair_atom1_atom2[j] = val2;
    }

    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j = 0; j < 5000; j++) {
      fscanf(fp1, "%le %le ", &val1, &val2);
      //(*eam_data)->radius[j]        = val1;
      (*eam_data)->pair_atom2_atom2[j] = val2;
    }

    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j = 0; j < 5000; j++) {
      fscanf(fp1, "%le %le ", &val1, &val2);
      //(*eam_data)->radius[j]        = val1;
      (*eam_data)->atom1_charge_Density[j] = val2;
    }

    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j = 0; j < 5000; j++) {
      fscanf(fp1, "%le %le ", &val1, &val2);
      //(*eam_data)->radius[j]        = val1;
      (*eam_data)->atom2_charge_Density[j]    = val2;
    }

    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j = 0; j < 5000; j++) {
      fscanf(fp1, "%le %le ", &val1, &val2);
      (*eam_data)->chargeDensity[j]           = val1;
      (*eam_data)->atom1_embedding_energy[j]  = val2;
    }

    fscanf(fp_init,"%s\n",str);
    fp1=fopen(str,"r");
    if(fp1==NULL)printf("flag Unable to open %s\n",str);

    for(j = 0; j < 5000; j++) {
      fscanf(fp1, "%le %le ", &val1, &val2);
      //(*eam_data)->radius[j]        = val1;
      (*eam_data)->atom2_embedding_energy[j] = val2;
    }
    (*eam_data)->min_eDen = (*eam_data)->chargeDensity[0];
    (*eam_data)->max_eDen = (*eam_data)->chargeDensity[4999];
}

struct radius_dependent_fields {
  int index;
  double radius;
  double p11;
  double p12;
  double p22;
  double eDen1;
  double eDen2;
}typedef rdf;

rdf* radius_retrive(binEAMpot* data, double radius) {
  if(radius < data->minRadius || radius > data->maxRadius){
    printf("Radius not in scope\n");
    return NULL;
  }

  int j = (int)((radius - data->minRadius)*5000)/(data->maxRadius - data->minRadius);
  rdf* new = (rdf*)malloc(sizeof(rdf));

  if(data->radius[j] < radius){
    while(data->radius[j] < radius)
      j++;
  }
  new->index = j;
  //interpolation functions
  new->p11   = linear_interpolator(data->radius[j-1], data->pair_atom1_atom1[j-1], data->radius[j], data->pair_atom1_atom1[j], radius);
  new->p12   = linear_interpolator(data->radius[j-1], data->pair_atom1_atom2[j-1], data->radius[j], data->pair_atom1_atom2[j], radius);
  new->p22   = linear_interpolator(data->radius[j-1], data->pair_atom2_atom2[j-1], data->radius[j], data->pair_atom2_atom2[j], radius);
  new->eDen1 = linear_interpolator(data->radius[j-1], data->atom1_charge_Density[j-1], data->radius[j], data->atom1_charge_Density[j], radius);
  new->eDen2 = linear_interpolator(data->radius[j-1], data->atom2_charge_Density[j-1], data->radius[j], data->atom2_charge_Density[j], radius);
  return new;
}

struct chargeDensity_dependent_fields {
  int index;
  double eDen;
  double embed1;
  double embed2;
}typedef eDen_df;

eDen_df* charge_density_retrive(binEAMpot* data, double eDen){
  if(eDen < data->min_eDen || eDen > data->max_eDen){
    printf("Radius not in scope\n");
    return NULL;
  }

  eDen_df* new;
  new = (eDen_df*)malloc(sizeof(eDen_df));
  int j = (int)((eDen - data->min_eDen)*5000)/(data->max_eDen - data->min_eDen);
  new->eDen   = eDen;
  new->embed1 = linear_interpolator(data->chargeDensity[j-1], data->atom1_embedding_energy[j-1], data->chargeDensity[j], data->atom1_embedding_energy[j], eDen);
  new->embed2 = linear_interpolator(data->chargeDensity[j-1], data->atom2_embedding_energy[j-1], data->chargeDensity[j], data->atom2_embedding_energy[j], eDen);

  return new;
}

void dataRetrivalTest() {
  double r = 4.123;
  binEAMpot* data = NULL;
  rdf*    nn = (rdf*)malloc(sizeof(rdf));
  eDen_df* a = (eDen_df*)malloc(sizeof(eDen_df));

  eam_data_read(&data, "file_list.txt", "Al", "Ni");
  nn = radius_retrive(data, r);
  printf("%d\n", nn->index);
  printf("%le\n", nn->p11);
  printf("%le\n", nn->p12);
  printf("%le\n", nn->p22);
  printf("%le\n", nn->eDen1);
  printf("%le\n", nn->eDen2);

  a = charge_density_retrive(data, 0.01);
  printf("\n");
  printf("%le\n", a->eDen);
  printf("%le\n", a->embed1);
  printf("%le\n", a->embed2);
}

void eam_data_read_test(){
  binEAMpot* data = NULL;
  eam_data_read(&data, "file_list.txt", "Al", "Ni");

  printf("%d\n", data->no_of_files);
  printf("%s\n", data->atom1);
  printf("%s\n", data->atom2);
}


int main() {
  //new_parameters_test();
  //eam_data_read_test();
  dataRetrivalTest();
}
