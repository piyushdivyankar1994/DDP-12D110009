/*!
   \file point3D.c
   \brief Source file for point3D.h
   \author Piyush Divyankar
   \date 01/09/2016
*/

#include "point3D.h"

/**
 * Allocates the memory, and creates a new point3D type.
 * @param  x x coordinate.
 * @param  y y coordinate
 * @param  z z coordinate
 * @return   pointer to newly created point3D type.
 */
point3D * point3D_newPoint(float x, float y, float z)
{
    point3D * new = (point3D *) malloc(sizeof(point3D));

    new->x = x;
    new->y = y;
    new->z = z;
    return new;
}

/**
 * Creates new point3D type with (0, 0, 0) as input
 * @return pointer to newly created point3D type.
 */
point3D * point3D_origin()
{
    point3D * new = (point3D *) malloc(sizeof(point3D));

    new->x = 0;
    new->y = 0;
    new->z = 0;
    return new;
}

/**
 * Adds two point3D types and creates a new point3D type to store the result
 * @param  v1 Pointer to first point3D to add.
 * @param  v2 Pointer to second point3D to add.
 * @return    Newly created point3D that holds the sum.
 */
point3D * point3D_addVectors(point3D * v1, point3D * v2)
{
    point3D * result = (point3D *) malloc(sizeof(point3D));

    result->x = v1->x + v2->x;
    result->y = v1->y + v2->y;
    result->z = v1->z + v2->z;
    return result;
}
/**
 * Compares two ::point3D objects and returns 1 if two are approximately equal
 * @param  v1 first point3D pointer
 * @param  v2 second point3D pointer
 * @return    1 if true -1 of false
 * FUTURE_CHANGES:10 Make it more intuitive like true false
 */
int point3D_isEqual(point3D * v1, point3D * v2)
{
    if ((v1->x - v2->x) < 1e-6 && (v1->y - v2->y) < 1e-6 && (v1->z - v2->z) < 1e-6)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

/**
 * Creates a new ::point3D object by subtacting first and second arguments.
 * @param  v1 first point3D pointer
 * @param  v2 second point3D pointer
 * @return    pointer to newly created point3D object.
 */
point3D * point3D_subtractVectors(point3D * v1, point3D * v2)
{
    point3D * result = (point3D *) malloc(sizeof(point3D));

    result->x = v1->x - v2->x;
    result->y = v1->y - v2->y;
    result->z = v1->z - v2->z;
    return result;
}
/**
 * Prints ::point3D data to STDOUT as "(x, y, z)".
 * @param a pointer to point3D object
 */
void point3D_dispPoint(point3D * a)
{
    if (a == NULL)
    {
        return;
    }
    printf("(%.2f\t%.2f\t%.2f)\t\n", a->x, a->y, a->z);
    return;
}
/**
 * Moves ::point3D by amount (\a +dx,\a +dy,\a +dz)
 * @param a  Pointer to point3D object
 * @param dx change in x-coordinate
 * @param dy change in y-coordinate
 * @param dz change in z-coordinate
 */
void point3DmoveTransform(point3D * a, float dx, float dy, float dz)
{
    a->x = a->x + dx;
    a->y = a->y + dx;
    a->z = a->z + dx;
}



/**
 * Converts a given array index in the range of Atomic crystal matrix and
 * parameters consistent with it it will return position location of atom at
 * that index in universal coordinate system of FCC lattice.
 * @param  index index from ::ATOM array
 * @param  p     parameter object pointer that is consistent with ATOM array
 * @return       pointer to newly created point3D object
 */
point3D * point3D_indexToPoint3D_fcc(int index, parameter * p)
{
    point3D * returnPoint = point3D_origin();
    int motiff_position = index % p->atoms_per_site;

    returnPoint->x = (float) ((index / p->atoms_per_site) % p->Nx);
    returnPoint->y = (float) ((index / (p->atoms_per_site * p->Nx)) % p->Ny);
    returnPoint->z = (float) ((index / (p->atoms_per_site * p->Nx * p->Ny)) % p->Nz);

    if (motiff_position == 0)
    {
        return returnPoint;
    }
    else if (motiff_position == 1)
    {
        point3DmoveTransform(returnPoint, 0.0, 0.5, 0.5);
        return returnPoint;
    }


    else if (motiff_position == 2)
    {
        point3DmoveTransform(returnPoint, 0.5, 0.0, 0.5);
        return returnPoint;
    }

    else
    {
        point3DmoveTransform(returnPoint, 0.5, 0.5, 0.0);
        return returnPoint;
    }
}

point3D * point3D_indexToPoint3D_bcc(int index, parameter *p) {
  point3D* returnPoint = point3D_origin();
  int k = index % 2;
  returnPoint->x = (float) ((index / p->atoms_per_site) % p->Nx);
  returnPoint->y = (float) ((index / (p->atoms_per_site * p->Nx)) % p->Ny);
  returnPoint->z = (float) ((index / (p->atoms_per_site * p->Nx * p->Ny)) % p->Nz);

  if(k == 1) {
    point3DmoveTransform(returnPoint, 0.5, 0.5, 0.5);
  }
  return returnPoint;
}
/**
 * Convertes a ::point3D data to a index in range of ::ATOM with respect to
 * ::parameter that is passed.
 * @param  a pointer to point3D object
 * @param  p pointer to parameter
 * @return   integer value that is the index of given point in array.
 * NOTE: Here it is assumed that crystal is FCC and there are certain assumptions
 * made about this crystal.
 * NOTE: Document how atom arrays work.
 */
int point3D_point3DtoIndexFCC(point3D * a, parameter * p)
{
    float fx = a->x - floor(a->x);
    float fy = a->y - floor(a->y);
    float fz = a->z - floor(a->z);

    int pos = ((int) floor(a->x) + p->Nx * (int) floor(a->y) + p->Nx * p->Ny * floor(a->z)) * 4;

    if (fx == 0 && fy == 0 && fz == 0)
    {
        return pos;
    }
    else if (fx == 0)
    {
        return pos + 1;
    }
    else if (fy == 0)
    {
        return pos + 2;
    }
    else if (fz == 0)
    {
        return pos + 3;
    }
    return pos;
}

int point3D_point3DtoIndexBCC(point3D * a, parameter *p)
{
    float fx = a->x - floor(a->x);
    float fy = a->y - floor(a->y);
    float fz = a->z - floor(a->z);

    int pos = ((int) floor(a->x) + p->Nx * (int) floor(a->y) + p->Nx * p->Ny * floor(a->z)) * 2;

    if (fx == 0 && fy == 0 && fz == 0)
    {
        return pos;
    }
    else {
        return pos+1;
    }
}
/**
 * For a given ::parameter and ::point3D applies the periodic boundary transform
 * such that it lies in the confined space of (\a Nx, \a Ny, \a Nz) as defined
 * in parameter.
 * @param k pointer to transform
 * @param p pointer to parameter
 */
void point3D_periodicBoundaryTransform(point3D * k, parameter * p)
{
    if (k->x < 0)
    {
        k->x = k->x + p->Nx;
    }
    else if (k->x > ((float) p->Nx - 0.5))
    {
        k->x = k->x - p->Nx;
    }
    if (k->y < 0)
    {
        k->y = k->y + p->Ny;
    }
    else if (k->y > ((float) p->Ny - 0.5))
    {
        k->y = k->y - p->Ny;
    }
    if (k->z < 0)
    {
        k->z = k->z + p->Nz;
    }
    else if (k->z > ((float) p->Nz - 0.5))
    {
        k->z = k->z - p->Nz;
    }
    // return k;
}

/**
 * Returns distance of ::point3D from origin in UC.
 * @param  a pointer to point3D
 * @return   distance from center or magnitude of the vector.
 */
float point3D_magnitude(point3D * a)
{
    return sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
}

/**
 * [point3D_distAtoB description]
 * @param  A [description]
 * @param  B [description]
 * @return   [description]
 */
float point3D_distAtoB(point3D * A, point3D * B)
{
    point3D * new = point3D_subtractVectors(A, B);
    float __retVal = point3D_magnitude(new);

    free(new);
    return __retVal;
}

/**
 * @breif Creates a index table that maps each index on the atomic lattice to its
 * neighbouring indices.
 * @details The atom at (0, 0, 0) in space is stored at index 0 of the ATOM array.
 * The neighbours of this are at (0.5, 0.5, 0) etc. These are stored in a
 * particular fashion and need to be transform from global coordinate system to
 * array index using the function point3D_point3DtoIndexFCC(..). This function is
 * costly and repetetive, i.e it does many repeat calculations.
 * @param  p     parameters of Simulation
 * @param  ngbrs coordinates of neighbouring fcc
 * @param  n     upto which nearest neighbours we want to compute
 * @return       pointer to a 2D array in form [siteNo][nearest_neighbour_no].
 *               the site number is which site we are looking at for neighbours,
 *               and nearest neighbours number is 0-11 for first nearest, 12-17
 *               for second nearest and so on.
 */
int* point3D_neighbourIndexTable_FCC(parameter *p, Sn_fcc *ngbrs, int n)
{
    int k = ngbrs->indices[n+1];
    int *ret_val = malloc(p->no_of_atoms * sizeof(int) * k);

    for(int i = 0; i < p->no_of_atoms; i++) {
        point3D * current = point3D_indexToPoint3D_fcc(i, p);
        for(int j = 0; j < k; j++) {
            point3D *nextNgbr = point3D_addVectors(current, &(ngbrs->s1n[j]));
            point3D_periodicBoundaryTransform(nextNgbr, p);
            int ngbrIndex = point3D_point3DtoIndexFCC(nextNgbr, p);
            ret_val[i*k + j] = ngbrIndex;
            free(nextNgbr);
        }
        free(current);
    }

    return ret_val;
}

int* point3D_neighbourIndexTable_BCC(parameter *p, Sn_bcc *ngbrs) {
    int * a = malloc(sizeof(int)*(p->no_of_atoms)*14);

    for(int i = 0; i < p->no_of_atoms; i++) {
        point3D * c = point3D_indexToPoint3D_bcc(i, p);
        for(int j = 0; j < 14; j++) {
            point3D * na = point3D_addVectors(c, &(ngbrs->s1n[j]));
            point3D_periodicBoundaryTransform(na, p);
            a[j + i * 14] = point3D_point3DtoIndexBCC(na, p);
            free(na);
        }
        free(c);
    }
    return a;
}
