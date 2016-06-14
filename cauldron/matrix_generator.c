#include<stdio.h>
#include<stdlib.h>
#include <gsl/gsl_rng.h>
#include<math.h>

#define N 10

double *parameter;

int read_parameter(char*,double*);

int main (void)
{
    const gsl_rng_type * T;
    gsl_rng * r;
    FILE *fp,*fw;
    fw=fopen("input_data.txt","w");
    char *ch;
    ch="parameters.txt";
    fp=fopen(ch,"r");
    int no_of_params;

    fscanf(fp,"%d ",&no_of_params);
    parameter=(double*)malloc(no_of_params*sizeof(double));

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    int n=read_parameter(ch,parameter);
    double u=parameter[0]* parameter[1]* parameter[2]; 
    
    double conc_g_p=parameter[5]/100;
    int lim_l=(int)parameter[2]*(100-conc_g_p)/200;
    int lim_u=(int)parameter[2]*(100+conc_g_p)/200;
    for(int i=0;i<u;i++)
    {
        double Nx=i%(int)parameter[0];
        double Ny=(i/(int)parameter[0])%(int)parameter[1];
        double Nz=(i/((int)parameter[0]*(int)parameter[1]))%(int)parameter[2];
       //NOTE Nickel Atoms reprented by 1
       //NOTE Aluminium atons represented by 0 
        if(Nz>=lim_l && Nz<=lim_u){
            fprintf(fw,"0 1 1 1\n");
        }
        else{
            for(int s=0;s<4;s++)
            {
                double k=gsl_rng_uniform(r);
                if(k<=conc_g_p)fprintf(fw,"1 ");
                else fprintf(fw,"0 ");

            }
            fprintf(fw,"\n");
        }
        
    }
    
    gsl_rng_free(r);
    fclose(fp);

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
       printf("%le\n",ret_val[i]);
    }
    return no_of_params;
}
               
