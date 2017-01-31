/*!
   \file chemicalPotential.c
   \brief Contains functions to calculate chemical potential.
   \author Piyush Divyankar
   \date 17/08/2016
 */

#include "chemicalPotential.h"

/*!
   \brief Return atom fraction of atom 1 in the atomic matrix
   \param a atomic matrix
   \pre Atomic matric must be supplied to the function and shouldn't be different
        in size than the supplied no_of_atoms.
   \return Atom fraction of atom 1
 */
float concentrationOfAtom1(int * a, int no_of_atoms)
{
    int count = 0;

    for (size_t i = 0; i < no_of_atoms; i++)
    {
        if (a[i] == 0)
        {
            count++;
        }
    }
    return (float) (count / no_of_atoms);
}

/*!
   \brief Returns no of atoms of type 1 in a atomic matrix
   \param a atomic matrix
   \param no_of_atoms total no of atoms in the matrix
   \param atom_key integer value that represents that atom in the matrix
   \return No of atoms of type 1 in atomic matrix
 */
int countOfAtom(int * a, int no_of_atoms, int atom_key)
{
    int count = 0;

    for (size_t i = 0; i < no_of_atoms; i++)
    {
        if (a[i] == atom_key)
        {
            count++;
        }
    }
    return count;
}

/*!
   \brief calculates chemical potential at a given index
   \param index
   \pre all relevant information i.e binary EAM potential, atomic matrix, parameter,
   neighbour information, should be read beforehand and supplied to this function
   \return Chemical potential at index.
 */
double chemicalPotentialAtIndex(int index, int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    double swap = energyToSwap(index, a, data, p, ngbrs);

    if (a[index] == 0)
    {
        double chemPot = swap * (-1) * (double) p->no_of_atoms;
        return chemPot;
    }
    else
    {
        double chemPot = swap * (double) p->no_of_atoms;
        return chemPot;
    }
}


double concentrationAtIndex(ATOM * atomicMatrix, int * ngbrTable, int index, int n)
{
    double count = 0;

    for (int i = 0; i < n; i++)
    {
        if (atomicMatrix[ngbrTable[index * n + i]] == 0)
        {
            count++;
        }
    }
    if (atomicMatrix[index] == 0)
    {
        count++;
    }
    return count / (double) n;
}

double * createConcentrationTable(ATOM * atomicMatrix, parameter * p, int * ngbrTable, int n)
{
    double * ret_val = malloc(sizeof(double) * p->no_of_atoms);

    for (int i = 0; i < p->no_of_atoms; i++)
    {
        ret_val[i] = concentrationAtIndex(atomicMatrix, ngbrTable, i, n);
    }
    return ret_val;
}

void printConcentrationTable(double * table, int initial_index, int final_index)
{
    for (int i = initial_index; i < final_index + 1; i++)
    {
        printf("%0.4f\n", table[i]);
    }
}

void updateConcentrationTable(double * table, int * ngbrIndex, int index, int n, int change)
{
    double dC = (change == 1 ? (-1/(double)(n+1)) : (1/(double)(n+1)));
    for(int i = 0; i < n; i++) {
        table[ngbrIndex[index * n + i]] += dC;
    }
    table[index] += dC;
}
