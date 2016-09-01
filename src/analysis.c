#include "analysis.h"

double analysis_totalEnergy(int *a, binEAMpot *data, parameter *p, Sn_fcc *ngbrs) {
  double ret_val = 0;

  for (size_t i = 0; i < p->no_of_atoms; i++) {
    ret_val = ret_val + energyAtIndexFCC(i, a, data, p, ngbrs);
  }

  return ret_val;
}

int orderedPhaseCount(int *a, parameter *p, Sn_fcc *ngbrs) {
  printf("oregnioegn\n");
  int count = 0;
  int j;
  for (size_t i = 0; i < p->no_of_atoms; i++) {
    point3D *new = point3D_indexToPoint3D_fcc(i, p);
    for (j = 0; j < 12; j++) {
      point3D *test = point3D_addVectors(new, &(ngbrs->s1n[j]));
      int t = point3D_point3DtoIndexFCC(test, p);
      if(a[t] == a[i]) {
        break;
      }
      free(test);
    }
    if(j == 12) {
      count++;
    }
    free(new);
  }
  printf("oregnioegn\n");
  return count;
}

int antiOrderedPhaseCount(int *a, parameter *p, Sn_fcc *fcc) {
  int count = 0;
  int j;
  int checkSite = 0;
  for (size_t i = 0; i < p->no_of_atoms; i++) {
    point3D *new = point3D_indexToPoint3D_fcc(i, p);
    checkSite = 0;
    for (j = 0; j < 12; j++) {
      point3D *test = point3D_addVectors(new, &(fcc->s1n[j]));
      int t = point3D_point3DtoIndexFCC(test, p);
      if(a[t] != a[i]) {
        checkSite++;
      }
    }
    if(checkSite == 11) {
      printf("flag %d \n",j);
      count++;
    }
  }
  return count;
}

void randomMatrixGeneratorFCC(parameter *p, char *ouput_file_name, unsigned long int seed) {
  const gsl_rng_type *T;
  gsl_rng *r;
  strcat(ouput_file_name, ".crystal.fcc");
  FILE *fp = fopen(ouput_file_name, "w");
  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (size_t i = 0; i < p->no_of_atoms; i++) {
    double u = gsl_rng_uniform(r);
    int a;
    if(u > 0.5)  a = 1;
    else a = 0;
    if(fwrite(&a, sizeof(int), 1, fp) == 1) continue;
    else i--;
  }
  fclose(fp);
  char str[1000];
  strcpy(str, ouput_file_name);
  strcat(str, ".random.gen.param");
  FILE *fr = fopen(str, "w");
  gsl_rng_fwrite(fr, r);
  //gsl_rng_free(r);
  //fclose(fr);
}

int totalAtomsInFile(FILE *fp) {
  fseek(fp, 0L , SEEK_END);
  int no_of_atoms = ftell(fp) / sizeof(int);
  fseek(fp, 0L, SEEK_SET);
  return no_of_atoms;
}

/// TODO:20 Binary type data storage that always talked about
