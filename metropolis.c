#include<stdio.h>
#include<stdlib.h>
#include <gsl/gsl_rng.h>
#include<math.h>


//This function takes difference in energy E and a scale factor that is by (default set to 1) and makes a random decission of the jump happening
//it returns 1 for making it happenin
//return 0 for not making it happen
/*int jump_decision(double,double);
      const gsl_rng_type * T;
      gsl_rng * r;

      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc (T); */


/*int jump_decision(double E,double scale_factor)
{
      
      const gsl_rng_type * T;
      gsl_rng * r;

      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);                                                          

      double u = gsl_rng_uniform (r);
      double probablity=exp((-1*E)/scale_factor);
      if(probablity>1)
          return 1;
      else if(probablity>u)
          return 1;
      else 
          return 0;

}*/

//driver function to test the code 

int main(void)
{ 
    double E=0;
    double scale_factor=1;
    FILE *fp;
    fp=fopen("testing_data_metropolis_code.dat","w");
    int i;
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);                                                          

    double u = gsl_rng_uniform (r);
    double probablity=exp((-1*E)/scale_factor);

    while(E<=1)
    {
        int count=0;
        double probablity=exp((-1*E)/scale_factor);
        printf("%le\n",probablity);
        for(i=0;i<1000;i++)
        {
            if(u<probablity)
                count++;
            u=gsl_rng_uniform(r);
        }
       
        fprintf(fp,"%le \t %le \t %le\n",E,((double)count/(double)i),probablity);
        E+=0.01;
    }
    gsl_rng_free(r); 
    return 0;
}




