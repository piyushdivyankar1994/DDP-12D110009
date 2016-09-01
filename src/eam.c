#include "eam.h"
#include "analysis.h"

int* atomicMatrixRead(char* fileName, parameter* p) {
	int* mat = (int*)malloc(sizeof(int) * p->no_of_atoms);
	FILE *fp = fopen(fileName, "r");
	if(fp == NULL) return NULL;
	for(int i = 0; i < p->no_of_atoms; i++){
		fscanf(fp, "%d", &mat[i]);
	}
	return mat;
}

// FUTURE_CHANGES:20 following function is alternate for atomicMatrixRead(...)

int* readCrystalFileFCC(char *fileName) {
  FILE *fp = fopen(fileName, "r");
  int n = totalAtomsInFile(fp);
  int *a = (int*)malloc(sizeof(int) * n);
  fread(a, sizeof(int), n, fp);
  return a;
}


void print_AtomicMatrix(int* mat, int begin, int end){
	for(int i = begin; i < end; i += 4) {
		printf("%d %d %d %d\n", mat[i], mat[i+1], mat[i+2], mat[i+3]);
	}
}

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

Sn_fcc* _defaultFCCNeighbours() {
	Sn_fcc* new = Sn_fcc_readNeighbours_fromFile("./neighbours/file_list_neighbours.txt");
	return new;
}

void eam_data_read(binEAMpot** eam_data, char *fileName, char atom1[2], char atom2[2])
{
    FILE *fp_init,*fp1;
    fp_init=fopen(fileName,"r");
    *eam_data = (binEAMpot*)malloc(sizeof(binEAMpot));
    //double* eam_data;
    //(*eam_dat)->atom1 = (char*)malloc(sizeof(atom1)+1);
    //(*eam_dat)->atom2 = (char*)malloc(sizeof(atom2)+1);

		///FIXME:0 assignment of the names of the atoms in EAM potential structure
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


rdf* rdf_radius_retrive(binEAMpot* data, double radius) {
  if(radius > data->maxRadius){
    rdf* new = (rdf*)malloc(sizeof(rdf));
		new->p11   = 0;
	  new->p12   = 0;
	  new->p22   = 0;
	  new->eDen1 = 0;
	  new->eDen2 = 0;
    return new;
  }

	else if (radius < data->minRadius) {
		rdf* new = (rdf*)malloc(sizeof(rdf));
		new->p11   = linear_interpolator(data->radius[0], data->pair_atom1_atom1[0], data->radius[1], data->pair_atom1_atom1[1], radius);
	  new->p12   = linear_interpolator(data->radius[0], data->pair_atom1_atom2[0], data->radius[1], data->pair_atom1_atom2[1], radius);
	  new->p22   = linear_interpolator(data->radius[0], data->pair_atom2_atom2[0], data->radius[1], data->pair_atom2_atom2[1], radius);
	  new->eDen1 = linear_interpolator(data->radius[0], data->atom1_charge_Density[0], data->radius[1], data->atom1_charge_Density[1], radius);
	  new->eDen2 = linear_interpolator(data->radius[0], data->atom2_charge_Density[0], data->radius[1], data->atom2_charge_Density[1], radius);
	}


  int j = (int)((radius - data->minRadius)*5000)/(data->maxRadius - data->minRadius);
  rdf* new = (rdf*)malloc(sizeof(rdf));

	if (j < 0) j = 0;
	else if (j > 5000) j = 5000;

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

eDen_df* eDen_df_charge_density_retrive(binEAMpot* data, double eDen){
	eDen_df *new = (eDen_df*)malloc(sizeof(eDen_df));
	if(eDen < data->min_eDen) {
		new->eDen   = eDen;
	  new->embed1 = linear_interpolator(data->chargeDensity[0], data->atom1_embedding_energy[0], data->chargeDensity[1], data->atom1_embedding_energy[1], eDen);
	  new->embed2 = linear_interpolator(data->chargeDensity[0], data->atom2_embedding_energy[0], data->chargeDensity[1], data->atom2_embedding_energy[1], eDen);
	}

	if(eDen > data->max_eDen) {
		new->eDen   = eDen;
	  new->embed1 = 0;
	  new->embed2 = 0;
		return new;
	}

	int j = (int)((eDen - data->min_eDen)*5000)/(data->max_eDen - data->min_eDen);

	if (j < 0) j = 0;
	else if (j > 5000) j = 5000;

  new->eDen   = eDen;
  new->embed1 = linear_interpolator(data->chargeDensity[j-1], data->atom1_embedding_energy[j-1], data->chargeDensity[j], data->atom1_embedding_energy[j], eDen);
  new->embed2 = linear_interpolator(data->chargeDensity[j-1], data->atom2_embedding_energy[j-1], data->chargeDensity[j], data->atom2_embedding_energy[j], eDen);

  return new;
}



double energyAtIndexFCC(int index, int* a, binEAMpot* data, parameter* p, Sn_fcc* ngbrs) {
	point3D* current = point3D_indexToPoint3D_fcc(index, p);
	double energy = 0;
	double chargeDen = 0;
	int noNgbrs = 134;

	for(int i = 0; i < noNgbrs; i++){
		point3D* k = point3D_addVectors(current, &(ngbrs->s1n[i]));
		double r = point3D_distAtoB(k, current);
		r = r * p->lattice_parameter;
		point3D_periodicBoundaryTransform(k, p);
		rdf* r_data = rdf_radius_retrive(data, r);
		int ngbrIndex = point3D_point3DtoIndexFCC(k, p);

		if (a[ngbrIndex] == a[index]) {
			if (a[ngbrIndex] == 0) {
				energy 		= energy + r_data->p11;
				chargeDen	= chargeDen + r_data->eDen1;
			}
			else {
				energy 		= energy + r_data->p22;
				chargeDen	= chargeDen + r_data->eDen2;
			}
		}
		else {
			if (a[ngbrIndex] == 0) {
				energy 		= energy + r_data->p12;
				chargeDen	= chargeDen + r_data->eDen1;
			}
			else {
				energy 		= energy + r_data->p12;
				chargeDen	= chargeDen + r_data->eDen2;
			}
		}
		free(r_data);
	}

	eDen_df* embeddingEnergy = eDen_df_charge_density_retrive(data, chargeDen);

	if (a[index] == 0) {
		energy = energy + embeddingEnergy->embed1;
	} else {
		energy = energy + embeddingEnergy->embed2;
	}
	free(embeddingEnergy);
	free(current);
	return energy;
}

void energyInMatrix(double **energyMatrix, int *a, binEAMpot *data, parameter *p, Sn_fcc *ngbrs) {
	 *energyMatrix = (double*)malloc(sizeof(double) * p->no_of_atoms);

	 for(int i = 0; i < p->no_of_atoms; i++) {
		 (*energyMatrix)[i] = energyAtIndexFCC(i, a, data, p, ngbrs);
	 }
}

void printEnergyMap(double *matrix, int init_index, int final_index) {
	for (int i = init_index; i < final_index; i++) {
		printf("#%d\t%le\n", i, matrix[i]);
	}
}

double energyToSwap(int index, int* a, binEAMpot* data, parameter* p, Sn_fcc* ngbrs) {
	point3D* current = point3D_indexToPoint3D_fcc(index, p);
	double energy1 = 0;
	double energy2 = 0;
	double chargeDen = 0;
	int noNgbrs = 134;

	for(int i = 0; i < noNgbrs; i++){
		point3D* k = point3D_addVectors(current, &(ngbrs->s1n[i]));
		double r = point3D_distAtoB(k, current);
		r = r * p->lattice_parameter;
		point3D_periodicBoundaryTransform(k, p);
		rdf* r_data = rdf_radius_retrive(data, r);
		int ngbrIndex = point3D_point3DtoIndexFCC(k, p);

		if (a[ngbrIndex] == a[index]) {
			if (a[ngbrIndex] == 0) {
				energy1 		= energy1 + r_data->p11;
				energy2			= energy2 + r_data->p12;
				chargeDen	= chargeDen + r_data->eDen1;
			}
			else {
				energy1 		= energy1 + r_data->p22;
				energy2			= energy2 + r_data->p12;
				chargeDen	= chargeDen + r_data->eDen2;
			}
		}
		else {
			if (a[ngbrIndex] == 0) {
				energy1 		= energy1 + r_data->p12;
				energy2			= energy2 + r_data->p22;
				chargeDen	= chargeDen + r_data->eDen1;
			}
			else {
				energy1 		= energy1 + r_data->p12;
				energy2			= energy2 + r_data->p11;
				chargeDen	= chargeDen + r_data->eDen2;
			}
		}
		free(r_data);
	}

	eDen_df* embeddingEnergy = eDen_df_charge_density_retrive(data, chargeDen);

	if (a[index] == 0) {
		energy1 = energy1 + embeddingEnergy->embed1;
		energy2 = energy2 + embeddingEnergy->embed2;
	} else {
		energy1 = energy1 + embeddingEnergy->embed2;
		energy2 = energy2 + embeddingEnergy->embed1;
	}
	free(embeddingEnergy);
	free(current);
	return energy2 - energy1;
}

void deltaEnergyMatrix(double **energyMatrix, int *a, binEAMpot *data, parameter *p, Sn_fcc *ngbrs) {
	 *energyMatrix = (double*)malloc(sizeof(double) * p->no_of_atoms);

	 for(int i = 0; i < p->no_of_atoms; i++) {
		 (*energyMatrix)[i] = energyToSwap(i, a, data, p, ngbrs);
	 }
}
