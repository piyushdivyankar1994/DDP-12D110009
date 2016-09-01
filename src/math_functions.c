// binary search function, searches the position of value that is lower limit of the smallest interval
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "math_functions.h"
// #include "eam_functions.h"
// #include "utils.h"

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

// exponentian function so we avoid using the default one ans save time
double exponential(double x)
{
    double sum = 1;

    for (int i = 25; i > 0; --i)
    {
        sum = 1 + x * sum / i;
    }
    return sum;
}

// Function to calculate euclidean distance between 2 points in 3D
double dist_3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

// Interpolates value of a variable given two points of a straight line
double linear_interpolator(double x1, double y1, double x2, double y2, double xc)
{
    double m = (y1 - y2) / (x1 - x2);   // slope

    return y1 + m * (xc - x1);
}
