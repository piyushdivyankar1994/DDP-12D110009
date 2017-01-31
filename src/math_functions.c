/*!
   \file math_functions.c
   \brief Source code for math_functions.h
   \author Piyush Divyankar
   \date 01/09/2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "math_functions.h"


/**
 * Binary search on a double array
 * @param  ptr pointer to array
 * @param  i   left index
 * @param  j   right index
 * @param  val value to search
 * @return     index of value just smaller than key value
 */
int binary_search(double * ptr, int i, int j, double val)
{
    if ((j - i) > 1)
    {
        int mid = (j + i) / 2;
        if (ptr[mid] == val)
        {
            return mid;
        }
        else if (ptr[mid] > val)
        {
            return binary_search(ptr, i, mid, val);
        }
        else
        {
            return binary_search(ptr, mid, j, val);
        }
    }
    else
    {
        return i;
    }
}

/**
 * Function for \f$e^x\f$
 * @param  x
 * @return   \f$e^x\f$
 */
double exponential(double x)
{
    double sum = 1;

    for (int i = 25; i > 0; --i)
    {
        sum = 1 + x * sum / i;
    }
    return sum;
}

/**
 * Distance between two points in 3D
 * @param  x1 [description]
 * @param  y1 [description]
 * @param  z1 [description]
 * @param  x2 [description]
 * @param  y2 [description]
 * @param  z2 [description]
 * @return    distane between \f$(x_1,y_1,z_1)\f$ and \f$(x_2,y_2,z_2)\f$
 * @deprecated This function replaced by ::point3D_distAtoB
 */
double dist_3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

/**
 * Interpolates linearly value between \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$
 * @param  x1
 * @param  y1
 * @param  x2
 * @param  y2
 * @param  xc
 * @return    return \f$y_c\f$ for \f$x_c\f$
 */
double linear_interpolator(double x1, double y1, double x2, double y2, double xc)
{
    double m = (y1 - y2) / (x1 - x2);   // slope

    return y1 + m * (xc - x1);
}
