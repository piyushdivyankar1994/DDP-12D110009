/*!
   \file math_functions.h
   \brief Function prototypes for general math functions.
   \author Piyush Divyankar
   \date 1/08/2016
 */

#ifndef _MATH_FUNCTIONS_H
  #define  _MATH_FUNCTIONS_H
  #include <stdio.h>
  #include <stdlib.h>
  #include <gsl/gsl_rng.h>
  #include <math.h>
  #include "math_functions.h"

int binary_search(double *, int, int, double);
double linear_interpolator(double, double, double, double, double);
double dist_3d(double, double, double, double, double, double);
double exponential(double);
#endif
