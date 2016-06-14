#include <stdio.h>
#include <stdlib.h>

#include<gsl/gsl_rng.h>

int main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;
  FILE *fp;

  fp=fopen("positions.dat","w");

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < n; i++) 
  {
      double u1 = gsl_rng_uniform (r)*3;
      double u2 = gsl_rng_uniform (r)*3;
	  fprintf(fp,"%lf %lf \n",u1,u2);
  }

  gsl_rng_free (r);

  return 0;
}
