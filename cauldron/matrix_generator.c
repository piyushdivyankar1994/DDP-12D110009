#include<stdio.h>
#include<stdlib.h>
#include <gsl/gsl_rng.h>
#include<math.h>
#include "eam_functions.h"
#include "math_functions.h"
#define N 10

//double *parameter;

int main (int argc, char **argv)
{
    argv++;
    char *data;
    data=*argv;
    random_crystal_generator("parameters.txt",data);
    return 0;
}

int read_parameter(char *filename,double ret_val[])
{
    FILE *fp;
    fp=fopen(filename,"r");
    int no_of_params;

    fscanf(fp,"%d ",&no_of_params);
    //ret_val=(double*)malloc(no_of_params*sizeof(double));
    for(int i=0;i<no_of_params;i++)
    {
        fscanf(fp,"%le ",&ret_val[i]);
    }
    return no_of_params;
}
