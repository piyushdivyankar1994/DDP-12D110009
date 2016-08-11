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
  int nearestNeighbours[7];
}typedef parameter;

/* moved to parameters.c
void print_parameters(parameter* a){
  printf("Size in x-direction           = %d\n", a->Nx);
  printf("Size in y-direction           = %d\n", a->Ny);
  printf("Size in z-direction           = %d\n", a->Nz);
  printf("Lattice Parameter             = %le\n", a->lattice_parameter);
  printf("No. of MonteCarlo Simulations = %d\n", a->N_MCS);
  printf("Temperature                   = %le\n", a->temperature);
  printf("Atoms per site                = %d\n", a->atoms_per_site);
  printf("Totol No of atoms             = %d\n", a->no_of_atoms );
}*/

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

void test_new_parameters(){
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

rdf* rdf_radius_retrive(binEAMpot* data, double radius) {
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

eDen_df* eDen_df_charge_density_retrive(binEAMpot* data, double eDen){
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

void test_dataRetrival() {
  double r = 4.123;
  binEAMpot* data = NULL;
  rdf*    nn = (rdf*)malloc(sizeof(rdf));
  eDen_df* a = (eDen_df*)malloc(sizeof(eDen_df));

  eam_data_read(&data, "file_list.txt", "Al", "Ni");
  nn = rdf_radius_retrive(data, r);
  printf("%d\n", nn->index);
  printf("%le\n", nn->p11);
  printf("%le\n", nn->p12);
  printf("%le\n", nn->p22);
  printf("%le\n", nn->eDen1);
  printf("%le\n", nn->eDen2);

  a = eDen_df_charge_density_retrive(data, 0.01);
  printf("\n");
  printf("%le\n", a->eDen);
  printf("%le\n", a->embed1);
  printf("%le\n", a->embed2);
}

void test_eam_data_read(){
  binEAMpot* data = NULL;
  eam_data_read(&data, "file_list.txt", "Al", "Ni");

  printf("%d\n", data->no_of_files);
  printf("%s\n", data->atom1);
  printf("%s\n", data->atom2);
}

struct point3D {
  float x;
  float y;
  float z;
}typedef point3D;

point3D* point3D_newPoint(float x, float y, float z){
	point3D* new = (point3D*)malloc(sizeof(point3D));
	new->x = x;
	new->y = y;
	new->z = z;
	return new;
}

point3D* point3D_origin(){
	point3D* new = (point3D*)malloc(sizeof(point3D));
	new->x = 0;
	new->y = 0;
	new->z = 0;
	return new;
}

point3D* point3D_addVectors(point3D* v1, point3D* v2){
	point3D* result = (point3D*)malloc(sizeof(point3D));
	result->x = v1->x + v2->x;
	result->y = v1->y + v2->y; 
	result->z = v1->z + v2->z;
	return result;
}

int point3D_isEqual(point3D* v1, point3D* v2) {
	if((v1->x - v2->x) < 1e-6 && (v1->y - v2->y) < 1e-6 && (v1->z - v2->z) < 1e-6) {
		return 1;
	}
	else return -1;
}

point3D* point3D_subtractVectors(point3D* v1, point3D* v2){
	point3D* result = (point3D*)malloc(sizeof(point3D));
	result->x = v1->x - v2->x;
	result->y = v1->y - v2->y; 
	result->z = v1->z - v2->z;
	return result;
}

void point3D_dispPoint(point3D* a){
	if(a == NULL) return;
	printf("(%.2f\t%.2f\t%.2f)\t\n", a->x, a->y, a->z);
	return;
}

point3D* point3D_indexToPoint3D_fcc(int index, parameter* p){
	point3D* returnPoint = point3D_origin();
	int motiff_position = index % p->atoms_per_site;          //type of atom
	
	returnPoint->x = (float)((index / p->atoms_per_site) % p->Nx);
	returnPoint->y = (float)((index / (p->atoms_per_site * p->Nx)) % p->Ny);
	returnPoint->z = (float)((index / (p->atoms_per_site * p->Nx * p->Ny)) % p->Nz);

	if(motiff_position == 0)                                        //if atom at first position(0,0,0) indicies are the global frame coordinates
	{
		return returnPoint;
	}
	else if (motiff_position == 1) {
		return point3D_addVectors(returnPoint, point3D_newPoint(0, 0.5, 0.5));
	}
	
	else if (motiff_position == 2) {
		return point3D_addVectors(returnPoint, point3D_newPoint(0.5, 0, 0.5));
	}
	
	else if (motiff_position == 3) {
		return point3D_addVectors(returnPoint, point3D_newPoint(0.5, 0.5, 0));
	}
	
	return returnPoint;
}

int point3D_point3DtoIndex(point3D* a, parameter* p){
	float fx = a->x - floor(a->x);
	float fy = a->y - floor(a->y);
	float fz = a->z - floor(a->z);
	
	int pos = ((int)floor(a->x) + p->Nx * (int)floor(a->y) + p->Nx * p->Ny * floor(a->z)) * 4;
	
	if(fx == 0 && fy == 0 && fz == 0){
		return pos;
	}
	else if(fx == 0)
		return pos + 1;
	else if(fy == 0)
		return pos + 2;
	else if(fz == 0)
		return pos + 3;
	return pos;
}

void point3D_periodicBoundaryTransform(point3D* k, parameter* p) {
	if(k->x < 0) {
		k->x = k->x + p->Nx;
	}
	else if(k->x > ((float)p->Nx - 0.5)){
		k->x = k->x - p->Nx;
	}
	if(k->y < 0) {
		k->y = k->y + p->Ny;
	}
	else if(k->y > ((float)p->Ny - 0.5)){
		k->y = k->y - p->Ny;
	}
	if(k->z < 0) {
		k->z = k->z + p->Nz;
	}
	else if(k->z > ((float)p->Nz - 0.5)){
		k->z = k->z - p->Nz;
	}
	//return k;
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

void test_point3D_indexToPoint3D_fcc(){
	parameter* p = _defaultFCCparameter();
	printf("Generating 10 random numbers from 0 to 3999 and displaying corresponding points\n");
	for(int i = 0; i < 15; i++){
		int index = rand() % p->no_of_atoms;
		printf("%d\t=\t", index);
		point3D_dispPoint(point3D_indexToPoint3D_fcc(index, p));
	}
}

float point3D_magnitude(point3D* a){
	return  sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
}

float point3D_distAtoB(point3D* A, point3D* B){
	return point3D_magnitude(point3D_subtractVectors(A, B));
}
void test_point3D(){
	point3D* a1 = point3D_origin();
	point3D_dispPoint(a1);
	a1 = point3D_newPoint(1, 1.5, 4.3);
	point3D_dispPoint(a1);
	point3D* a2 = point3D_newPoint(1.1, 3.4, 5.3);
	point3D_dispPoint(a2);
	point3D* sum = point3D_addVectors(a1, a2);
	point3D_dispPoint(sum);
	point3D* diff = point3D_subtractVectors(a1, a2);
	point3D_dispPoint(diff);
	printf("Distance Test from point3D_origin\n");
	a1 = point3D_newPoint(1, 1, 1);
	printf("%f\n", point3D_magnitude(a1));
	a1 = point3D_newPoint(1, 1.5, 1);
	printf("%f\n", point3D_magnitude(a1));
	a1 = point3D_newPoint(1.4, 1, 1);
	printf("%f\n", point3D_magnitude(a1));
	a1 = point3D_newPoint(0, 0, 5);
	printf("%f\n", point3D_magnitude(a1));
	
	printf("Distance between two points\n");
	for(int i = 0; i < 3; i++){
		a1 = point3D_newPoint(rand() / (float) RAND_MAX, rand() / (float) RAND_MAX, rand() / (float) RAND_MAX);
		a2 = point3D_newPoint(rand() / (float) RAND_MAX, rand() / (float) RAND_MAX, rand() / (float) RAND_MAX);
		point3D_dispPoint(a1);
		point3D_dispPoint(a2);
		printf("Distance = %f\n\n", point3D_distAtoB(a1, a2));
		
	}
	
	return;
}

void test_defaultFCCparameter(){
	parameter* new = _defaultFCCparameter();
	print_parameters(new);	
}

struct neighbours_fcc{
	point3D s1n[12];
	point3D s2n[6];
	point3D s3n[24];
	point3D s4n[12];
	point3D s5n[24];
	point3D s6n[8];
	point3D s7n[48];
	int indices[7];
}typedef Sn_fcc;

Sn_fcc* Sn_fcc_readNeighbours_fromFile(char* fileList){
	FILE* fp;
	Sn_fcc* new = (Sn_fcc*)malloc(sizeof(Sn_fcc));
	new->indices[0] = 12;
	new->indices[1] = 6;
	new->indices[2] = 24;
	new->indices[3] = 12;
	new->indices[4] = 24;
	new->indices[5] = 8;
	new->indices[6] = 48;
	fp = fopen(fileList, "r");
	if(fp == NULL) printf("flag\n");
	char str[35];
	for(int i = 0; i < 7; i++){
		fscanf(fp, "%s", str);
		FILE* fr = fopen(str , "r");
		int k = new->indices[i];
		
		if(fr == NULL)printf("flag\n");
		
		for(int j = 0; j < k; j++){
			float x, y, z;
			fscanf(fr, "%f", &x);
			fscanf(fr, "%f", &y);
			fscanf(fr, "%f", &z);
			
			if(i == 0) {
				new->s1n[j].x = x;
				new->s1n[j].y = y;
				new->s1n[j].z = z;
				}
			if(i == 1) {
				new->s2n[j].x = x;
				new->s2n[j].y = y;
				new->s2n[j].z = z;
				}
			if(i == 2) {
				new->s3n[j].x = x;
				new->s3n[j].y = y;
				new->s3n[j].z = z;
				}
			if(i == 3) {
				new->s4n[j].x = x;
				new->s4n[j].y = y;
				new->s4n[j].z = z;
				}
			if(i == 4) {
				new->s5n[j].x = x;
				new->s5n[j].y = y;
				new->s5n[j].z = z;
				}
			if(i == 5) {
				new->s6n[j].x = x;
				new->s6n[j].y = y;
				new->s6n[j].z = z;
				}
			if(i == 6) {
				new->s7n[j].x = x;
				new->s7n[j].y = y;
				new->s7n[j].z = z;
				}
		}
	}
	return new;
}

void print_Neighbours(Sn_fcc* a){
	for(int i = 0; i < 7; i++){
		int k = a->indices[i];
		for(int j = 0; j < k; j++){
			switch(i){
				case 0:
					point3D_dispPoint(&(a->s1n[j]));
					break;
				case 1:
					point3D_dispPoint(&(a->s2n[j]));
					break;
				case 2:
					point3D_dispPoint(&(a->s3n[j]));
					break;
				case 3:
					point3D_dispPoint(&(a->s4n[j]));
					break;
				case 4:
					point3D_dispPoint(&(a->s5n[j]));
					break;
				case 5:
					point3D_dispPoint(&(a->s6n[j]));
					break;
				case 6:
					point3D_dispPoint(&(a->s7n[j]));
					break;
				
			}
		}
	}
}

void test_Sn_fcc_readNeighbours_fromFile(){
	Sn_fcc* new = Sn_fcc_readNeighbours_fromFile("file_list_neighbours.txt");
	//print_Neighbours(new);
	printf("Everything is continuoulsy stored so if i goes from 0 \
	to sum(indicies) then s1n to s7n can be accessed as s1n[i] \
	this is demonstrated by printing 13 points\n ");
	for(int i = 0; i < 13; i++){
		point3D_dispPoint(&(new->s1n[i]));
	}
}

Sn_fcc* _defaultFCCNeighbours() {
	Sn_fcc* new = Sn_fcc_readNeighbours_fromFile("file_list_neighbours.txt");
	return new;
}

int* atomicMatrixRead(char* fileName, parameter* p) {
	int* mat = (int*)malloc(sizeof(int) * p->no_of_atoms);
	FILE *fp = fopen(fileName, "r");
	if(fp == NULL) return NULL;
	for(int i = 0; i < p->no_of_atoms; i++){
		fscanf(fp, "%d", &mat[i]);
	}
	return mat;
}

void print_AtomicMatrix(int* mat, int begin, int end){
	for(int i = begin; i < end; i += 4) {
		printf("%d %d %d %d\n", mat[i], mat[i+1], mat[i+2], mat[i+3]);
	}
}


void test_AtomicMatrixRead(){
	parameter* p = _defaultFCCparameter();
	int *a = atomicMatrixRead("out.txt", p);
	print_AtomicMatrix(a, 20, 50);
}

int hash_Radius(double value, binEAMpot* data){
}

int hash_eDen(double value, binEAMpot* data){
}

double energyAtIndexFCC(int index, int* a, binEAMpot* data, parameter* p, Sn_fcc* ngbrs) {
	point3D* current = point3D_indexToPoint3D_fcc(index, p);
	point3D_dispPoint(current);
	double energy = 0;
	double chargeDen = 0;
	int noNgbrs = 134;
	
	for(int i = 0; i < noNgbrs; i++){
		point3D* k = point3D_addVectors(current, &(ngbrs->s1n[i]));
		double r = point3D_distAtoB(k, current);
		//k = point3D_periodicBoundaryTransform(k, p);
		rdf* r_data = rdf_radius_retrive(data, r);
		int ngbrIndex = point3D_point3DtoIndex(k, p);
		if(a[index] == 0){
			if(a[ngbrIndex] == 0){
				
			}
		}
	}
	
	return energy;
}

void test_energyAtIndexFCC() { 
	int index = 100;
	binEAMpot* data = NULL;
	eam_data_read(&data, "file_list.txt", "Al", "Ni");
	parameter* p = _defaultFCCparameter();
	int *a = atomicMatrixRead("out.txt", p);
	//test_AtomicMatrixRead();
	Sn_fcc* new = _defaultFCCNeighbours();
	double e = energyAtIndexFCC(index, a, data, p, new);
	printf("at index = %d energy = %f", index, e);
	
}

void test_point3D_point3DtoIndex(){
	parameter* p = _defaultFCCparameter();
	point3D* p1 = point3D_newPoint(2, 0, 0);
	point3D* p2 = point3D_newPoint(0.5, 0.5, 1);
	point3D* p3 = point3D_newPoint(2.0, 0.5, 0.5);
	point3D* p4 = point3D_newPoint(0.5, 1, 0.5);
	
	printf("%d\t", point3D_point3DtoIndex(p1, p));
	point3D_dispPoint(p1);
	printf("%d\t", point3D_point3DtoIndex(p2, p));
	point3D_dispPoint(p2);
	printf("%d\t", point3D_point3DtoIndex(p3, p));
	point3D_dispPoint(p3);
	printf("%d\t", point3D_point3DtoIndex(p4, p));
	point3D_dispPoint(p4);
}

void test_point3D_periodicBoundaryTransform(){
	point3D* new = point3D_newPoint(-10.5, 11.5, 16.5);					
	parameter* p = _defaultFCCparameter();
	printf("test function for periodic boundary transform\n");
	point3D_dispPoint(new);
	//point3D* test = 
	point3D_periodicBoundaryTransform(new, p);
	point3D_dispPoint(new);
}

int main() {
  //test_new_parameters();
  //test_eam_data_read();
  //test_dataRetrival();
  //test_point3D();
  //test_defaultFCCparameter();
  //test_point3D_indexToPoint3D_fcc();
  //test_Sn_fcc_readNeighbours_fromFile();
  //test_AtomicMatrixRead();
  //test_energyAtIndex();
  //test_point3D_point3DtoIndex();
  test_point3D_periodicBoundaryTransform();
}
