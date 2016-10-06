/*!
   \file eam.c
   \brief Source file for ean.h
   \author Piyush Divyankar
   \date 01/09/2016
 */

#include "eam.h"
#include "analysis.h"

/**
 * Reads a atomic matrix in plain text format
 * @param  fileName FILE of atomic matrix in plain text
 * @param  p        ::parameter for simulation
 * @return          pointer to created atomic matrix
 */
int * atomicMatrixRead(char * fileName, parameter * p)
{
    ATOM * mat = (int *) malloc(sizeof(int) * p->no_of_atoms);
    FILE * fp = fopen(fileName, "r");

    if (fp == NULL)
    {
        return NULL;
    }
    for (int i = 0; i < p->no_of_atoms; i++)
    {
        fscanf(fp, "%d", &mat[i]);
    }
    return mat;
}

// DONE:20 following function is alternate for atomicMatrixRead(...)
/**
 * Reads a atomic matrix that was saved using .crystal.fcc format
 * @param  fileName Name of th file
 * @return          Pointer to ::ATOM array
 */
int * readCrystalFileFCC(char * fileName)
{
    FILE * fp = fopen(fileName, "r");
    int n = totalAtomsInFile(fp);
    ATOM * a = (int *) malloc(sizeof(int) * n);

    fread(a, sizeof(int), n, fp);
    return a;
}

/**
 * Prints atomic matrix atoms to STDOUT in a range
 * @param mat   ::ATOM matrix
 * @param begin start printing at this index
 * @param end   stop at this index
 */
void print_AtomicMatrix(ATOM * mat, int begin, int end)
{
    for (int i = begin; i < end; i += 4)
    {
        printf("%d %d %d %d\n", mat[i], mat[i + 1], mat[i + 2], mat[i + 3]);
    }
}

// TODO Write more documentation for this
// TODO Very bad code needs to be more concise.
/**
 * Reads ::neighbours_fcc from a file that contains list of files that conatin
 * list of neighbours.
 * @param  fileList name of file containing relevant lsit
 * @return          pointer to neighbours_fcc structure that is created
 */
Sn_fcc * Sn_fcc_readNeighbours_fromFile(char * fileList)
{
    FILE * fp;
    Sn_fcc * new = (Sn_fcc *) malloc(sizeof(Sn_fcc));

    new->indices[0] = 12;
    new->indices[1] = 6;
    new->indices[2] = 24;
    new->indices[3] = 12;
    new->indices[4] = 24;
    new->indices[5] = 8;
    new->indices[6] = 48;
    fp = fopen(fileList, "r");
    if (fp == NULL)
    {
        printf("flag\n");
    }
    char str[35];
    for (int i = 0; i < 7; i++)
    {
        fscanf(fp, "%s", str);
        FILE * fr = fopen(str, "r");
        int k = new->indices[i];

        if (fr == NULL)
        {
            printf("flag\n");
        }

        for (int j = 0; j < k; j++)
        {
            float x, y, z;
            fscanf(fr, "%f", &x);
            fscanf(fr, "%f", &y);
            fscanf(fr, "%f", &z);

            if (i == 0)
            {
                new->s1n[j].x = x;
                new->s1n[j].y = y;
                new->s1n[j].z = z;
            }
            if (i == 1)
            {
                new->s2n[j].x = x;
                new->s2n[j].y = y;
                new->s2n[j].z = z;
            }
            if (i == 2)
            {
                new->s3n[j].x = x;
                new->s3n[j].y = y;
                new->s3n[j].z = z;
            }
            if (i == 3)
            {
                new->s4n[j].x = x;
                new->s4n[j].y = y;
                new->s4n[j].z = z;
            }
            if (i == 4)
            {
                new->s5n[j].x = x;
                new->s5n[j].y = y;
                new->s5n[j].z = z;
            }
            if (i == 5)
            {
                new->s6n[j].x = x;
                new->s6n[j].y = y;
                new->s6n[j].z = z;
            }
            if (i == 6)
            {
                new->s7n[j].x = x;
                new->s7n[j].y = y;
                new->s7n[j].z = z;
            }
        }
    }
    return new;
}

/**
 * Prints neighbours to STDOUT as a list of ...(x, y, z) ... points
 * @param a pointer to neighbour structure
 */
void print_Neighbours(Sn_fcc * a)
{
    for (int i = 0; i < 7; i++)
    {
        int k = a->indices[i];
        for (int j = 0; j < k; j++)
        {
            switch (i)
            {
                case 0:
                    point3D_dispPoint(&(a->s1n[j]));
                    break;
                case 1:
                    point3D_dispPoint(&(a->s2n[j]));
                    break;
                case 2:
                    point3D_dispPoint(&(a->s3n[j]));
                    break;
                case 3:
                    point3D_dispPoint(&(a->s4n[j]));
                    break;
                case 4:
                    point3D_dispPoint(&(a->s5n[j]));
                    break;
                case 5:
                    point3D_dispPoint(&(a->s6n[j]));
                    break;
                case 6:
                    point3D_dispPoint(&(a->s7n[j]));
                    break;

            }
        }
    }
}

/**
 * Sets neighbours to default as FCC crystal
 * @return pointer to type
 */
Sn_fcc * _defaultFCCNeighbours()
{
    Sn_fcc * new = Sn_fcc_readNeighbours_fromFile("./neighbours/file_list_neighbours.txt");

    return new;
}

/**
 * Reads EAM data from a fixed file EAM_Ni_Al/file_list.txt
 * @param fileName name of the file
 * @param atom1    IUPAC symbol for atom 1.
 * @param atom2    IUPAC symbol for atom 2.
 */
binEAMpot * eam_data_read(char * fileName, char atom1[2], char atom2[2])
{
    FILE * fp_init, * fp1;

    fp_init = fopen(fileName, "r");
    binEAMpot * eam_data = (binEAMpot *) malloc(sizeof(binEAMpot));
    // double* eam_data;
    // (*eam_dat)->atom1 = (char*)malloc(sizeof(atom1)+1);
    // (*eam_dat)->atom2 = (char*)malloc(sizeof(atom2)+1);

    // /FIXME:10 assignment of the names of the atoms in EAM potential structure
    (eam_data)->atom1[0] = atom1[0];
    (eam_data)->atom1[1] = atom1[1];

    (eam_data)->atom2[0] = atom2[0];
    (eam_data)->atom2[1] = atom2[1];
    (eam_data)->no_of_files = 7;

    char str[100];
    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    double val1, val2;
    int j;
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        (eam_data)->radius[j]        = val1;
        (eam_data)->pair_atom1_atom1[j] = val2;
    }
    (eam_data)->minRadius = (eam_data)->radius[0];
    (eam_data)->maxRadius = (eam_data)->radius[4999];

    // printf("%d %le\n", j, (*eam_data)->radius[234]);
    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        // (*eam_data)->radius[j]        = val1;
        (eam_data)->pair_atom1_atom2[j] = val2;
    }

    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        // (*eam_data)->radius[j]        = val1;
        (eam_data)->pair_atom2_atom2[j] = val2;
    }

    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        // (*eam_data)->radius[j]        = val1;
        (eam_data)->atom1_charge_Density[j] = val2;
    }

    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        // (*eam_data)->radius[j]        = val1;
        (eam_data)->atom2_charge_Density[j]    = val2;
    }

    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        (eam_data)->chargeDensity[j]           = val1;
        (eam_data)->atom1_embedding_energy[j]  = val2;
    }

    fscanf(fp_init, "%s\n", str);
    fp1 = fopen(str, "r");
    if (fp1 == NULL)
    {
        printf("flag Unable to open %s\n", str);
    }

    for (j = 0; j < 5000; j++)
    {
        fscanf(fp1, "%le %le ", &val1, &val2);
        // (*eam_data)->radius[j]        = val1;
        (eam_data)->atom2_embedding_energy[j] = val2;
    }
    (eam_data)->min_eDen = (eam_data)->chargeDensity[0];
    (eam_data)->max_eDen = (eam_data)->chargeDensity[4999];
    return eam_data;
}

/**
 * In the EAM tables we have radius values common in 5 quantities as mentioned
 * in ::radius_dependent_fields. This function takes a key value and using a hash
 * function determines the index of that value in table, and linearly interpolates
 * other radius dependent quantities
 * @param  data   pointer to EAM table
 * @param  radius key value
 * @return        pointer to newly created ::rdf
 */
rdf * rdf_radius_retrive(binEAMpot * data, double radius)
{

    rdf * new = (rdf *) malloc(sizeof(rdf));

    new->radius = radius;
    if (radius > data->maxRadius)
    {
        new->p11   = 0;
        new->p12   = 0;
        new->p22   = 0;
        new->eDen1 = 0;
        new->eDen2 = 0;
        return new;
    }

    else if (radius < data->minRadius)
    {
        new->p11   = linear_interpolator(data->radius[0], data->pair_atom1_atom1[0], data->radius[1], data->pair_atom1_atom1[1], radius);
        new->p12   = linear_interpolator(data->radius[0], data->pair_atom1_atom2[0], data->radius[1], data->pair_atom1_atom2[1], radius);
        new->p22   = linear_interpolator(data->radius[0], data->pair_atom2_atom2[0], data->radius[1], data->pair_atom2_atom2[1], radius);
        new->eDen1 = linear_interpolator(data->radius[0], data->atom1_charge_Density[0], data->radius[1], data->atom1_charge_Density[1], radius);
        new->eDen2 = linear_interpolator(data->radius[0], data->atom2_charge_Density[0], data->radius[1], data->atom2_charge_Density[1], radius);
    }


    int j = (int) ((radius - data->minRadius) * 5000) / (data->maxRadius - data->minRadius);


    if (j < 0)
    {
        j = 0;
    }
    else if (j > 5000)
    {
        j = 5000;
    }

    if (data->radius[j] < radius)
    {
        while (data->radius[j] < radius)
        {
            j++;
        }
    }
    new->index = j;
    // interpolation functions
    new->p11   = linear_interpolator(data->radius[j - 1], data->pair_atom1_atom1[j - 1], data->radius[j], data->pair_atom1_atom1[j], radius);
    new->p12   = linear_interpolator(data->radius[j - 1], data->pair_atom1_atom2[j - 1], data->radius[j], data->pair_atom1_atom2[j], radius);
    new->p22   = linear_interpolator(data->radius[j - 1], data->pair_atom2_atom2[j - 1], data->radius[j], data->pair_atom2_atom2[j], radius);
    new->eDen1 = linear_interpolator(data->radius[j - 1], data->atom1_charge_Density[j - 1], data->radius[j], data->atom1_charge_Density[j], radius);
    new->eDen2 = linear_interpolator(data->radius[j - 1], data->atom2_charge_Density[j - 1], data->radius[j], data->atom2_charge_Density[j], radius);
    return new;
}
/**
 * In the EAM tables we have embedding energy that depends on charge density.
 * This function takes a key value and using a hash function determines the
 * index of that value in table, and linearly interpolates
 * other charge density dependent quantities.
 * @param  data   pointer to EAM table
 * @param  radius key value
 * @return        pointer to newly created ::eDen_df
 */
eDen_df * eDen_df_charge_density_retrive(binEAMpot * data, double eDen)
{
    eDen_df * new = (eDen_df *) malloc(sizeof(eDen_df));

    if (eDen < data->min_eDen)
    {
        new->eDen   = eDen;
        new->embed1 = linear_interpolator(data->chargeDensity[0], data->atom1_embedding_energy[0], data->chargeDensity[1], data->atom1_embedding_energy[1], eDen);
        new->embed2 = linear_interpolator(data->chargeDensity[0], data->atom2_embedding_energy[0], data->chargeDensity[1], data->atom2_embedding_energy[1], eDen);
    }

    if (eDen > data->max_eDen)
    {
        new->eDen   = eDen;
        new->embed1 = 0;
        new->embed2 = 0;
        return new;
    }

    int j = (int) ((eDen - data->min_eDen) * 5000) / (data->max_eDen - data->min_eDen);

    if (j < 0)
    {
        j = 0;
    }
    else if (j > 5000)
    {
        j = 5000;
    }

    new->eDen   = eDen;
    new->embed1 = linear_interpolator(data->chargeDensity[j - 1], data->atom1_embedding_energy[j - 1], data->chargeDensity[j], data->atom1_embedding_energy[j], eDen);
    new->embed2 = linear_interpolator(data->chargeDensity[j - 1], data->atom2_embedding_energy[j - 1], data->chargeDensity[j], data->atom2_embedding_energy[j], eDen);

    return new;
}

/**
 * Calculates energy at a given array index in ATOM array.
 * @param  index evaluated at this index.
 * @param  a     Atomic matrix
 * @param  data  EAM potential data.
 * @param  p     parameters
 * @param  ngbrs List of neighbours.
 * @return       value of energy
 * @note         This function has been tested in ::test_energyAtIndexFCC. Refer it
 *               for usage information.
 */
double energyAtIndexFCC(int index, int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    point3D * current = point3D_indexToPoint3D_fcc(index, p);
    double energy = 0;
    double chargeDen = 0;
    int noNgbrs = 134;

    for (int i = 0; i < noNgbrs; i++)
    {
        point3D * k = point3D_addVectors(current, &(ngbrs->s1n[i]));
        double r = point3D_distAtoB(k, current);
        r = r * p->lattice_parameter;
        point3D_periodicBoundaryTransform(k, p);
        rdf * r_data = rdf_radius_retrive(data, r);
        int ngbrIndex = point3D_point3DtoIndexFCC(k, p);

        if (a[ngbrIndex] == a[index])
        {
            if (a[ngbrIndex] == 0)
            {
                energy          = energy + r_data->p11;
                chargeDen       = chargeDen + r_data->eDen1;
            }
            else
            {
                energy          = energy + r_data->p22;
                chargeDen       = chargeDen + r_data->eDen2;
            }
        }
        else
        {
            if (a[ngbrIndex] == 0)
            {
                energy          = energy + r_data->p12;
                chargeDen       = chargeDen + r_data->eDen1;
            }
            else
            {
                energy          = energy + r_data->p12;
                chargeDen       = chargeDen + r_data->eDen2;
            }
        }
        free(r_data);
        free(k);
    }

    eDen_df * embeddingEnergy = eDen_df_charge_density_retrive(data, chargeDen);

    if (a[index] == 0)
    {
        energy = energy + embeddingEnergy->embed1;
    }
    else
    {
        energy = energy + embeddingEnergy->embed2;
    }
    free(embeddingEnergy);
    free(current);
    return energy;
}

double energyAtIndexFCC_ver2(int index, int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs, int noNgbrs)
{
    point3D * current = point3D_indexToPoint3D_fcc(index, p);
    double energy = 0;
    double chargeDen = 0;

    for (int i = 0; i < noNgbrs; i++)
    {
        point3D * k = point3D_addVectors(current, &(ngbrs->s1n[i]));
        double r = point3D_distAtoB(k, current);
        r = r * p->lattice_parameter;
        point3D_periodicBoundaryTransform(k, p);
        rdf * r_data = rdf_radius_retrive(data, r);
        int ngbrIndex = point3D_point3DtoIndexFCC(k, p);

        if (a[ngbrIndex] == a[index])
        {
            if (a[ngbrIndex] == 0)
            {
                energy          = energy + r_data->p11;
                chargeDen       = chargeDen + r_data->eDen1;
            }
            else
            {
                energy          = energy + r_data->p22;
                chargeDen       = chargeDen + r_data->eDen2;
            }
        }
        else
        {
            if (a[ngbrIndex] == 0)
            {
                energy          = energy + r_data->p12;
                chargeDen       = chargeDen + r_data->eDen1;
            }
            else
            {
                energy          = energy + r_data->p12;
                chargeDen       = chargeDen + r_data->eDen2;
            }
        }
        free(r_data);
        free(k);
    }

    eDen_df * embeddingEnergy = eDen_df_charge_density_retrive(data, chargeDen);

    if (a[index] == 0)
    {
        energy = energy + embeddingEnergy->embed1;
    }
    else
    {
        energy = energy + embeddingEnergy->embed2;
    }
    free(embeddingEnergy);
    free(current);
    return energy;
}

/**
 * Looks up the energy for atom at a given array index in ATOM array, from
 * energyTableInstantLookup[.][.][.][.]
 * @param  index evaluated at this index.
 * @param  a     Atomic matrix
 * @param  data  EAM potential data.
 * @param  p     parameters
 * @param  ngbrs List of neighbours.
 * @return       value of energy
 * @note         This function has been tested in ::test_energyAtIndexFCC_fast.
 *               Refer it for usage information.
 */

double energyAtIndexFCC_fast(int index, int * a, parameter * p, Sn_fcc * ngbrs)
{
    point3D * current = point3D_indexToPoint3D_fcc(index, p);
    int Ni[3] = {0};
    for (size_t n = 0; n < 3; n++)
    {
        for (int i = 0; i < ngbrs->indices[n]; i++)
        {
            point3D * k = point3D_addVectors(current, &(ngbrs->s1n[i]));

            point3D_periodicBoundaryTransform(k, p);
            int ngbrIndex = point3D_point3DtoIndexFCC(k, p);
            if (a[ngbrIndex] == 1)
            {
                Ni[n]++;
            }
            free(k);
        }
    }
    if (a[index] == 0)
    {
        return energyTableInstantLookup[0][Ni[0]][Ni[1]][Ni[2]];
    }
    else
    {
        return energyTableInstantLookup[1][Ni[0]][Ni[1]][Ni[2]];
    }
}



/**
 * Calculates bond energy of all atoms in a given ATOM array, and stores it in a
 * double array, at corresponding index with that of ATOM array.
 * @param  a     Atomic matrix
 * @param  data  EAM potential data.
 * @param  p     parameters
 * @param  ngbrs List of neighbours.
 * @return       value of energy
 * @note         This function has been tested in ::test_energyAtIndexFCC. Refer it
 *               for usage information.
 */
double * energyInMatrix(int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    double * energyMatrix = (double *) malloc(sizeof(double) * p->no_of_atoms);

    for (int i = 0; i < p->no_of_atoms; i++)
    {
        (energyMatrix)[i] = energyAtIndexFCC(i, a, data, p, ngbrs);
    }
    return energyMatrix;
}

double * energyInMatrix_ver2(int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs, int noNgbrs)
{
    double * energyMatrix = (double *) malloc(sizeof(double) * p->no_of_atoms);

    for (int i = 0; i < p->no_of_atoms; i++)
    {
        (energyMatrix)[i] = energyAtIndexFCC_ver2(i, a, data, p, ngbrs, noNgbrs);
    }
    return energyMatrix;
}

/**
 * Prints a double array to STDOUT.
 * @param matrix      double array
 * @param init_index  begin print at this index
 * @param final_index end at this index
 * @note              Created for use in context of energy matrix.
 */
void printEnergyMap(double * matrix, int init_index, int final_index)
{
    for (int i = init_index; i < final_index; i++)
    {
        printf("#%d\t%le\n", i, matrix[i]);
    }
}

/**
 * Function to calculate the energy required to swap an atom at a given index.
 * If say that the atom was initally atom1 then on swapping it would be atom2.
 *
 * @param  index index of ATOM array
 * @param  a     ATOM array
 * @param  data  EAM potential
 * @param  p     parameters of the system
 * @param  ngbrs list of neighbour atoms.
 * @return       change in energy when atom at the given index is swapped
 */
double energyToSwap(int index, int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    point3D * current = point3D_indexToPoint3D_fcc(index, p);
    double energy1 = 0;
    double energy2 = 0;
    double chargeDen = 0;
    int noNgbrs = 134;

    for (int i = 0; i < noNgbrs; i++)
    {
        point3D * k = point3D_addVectors(current, &(ngbrs->s1n[i]));
        double r = point3D_distAtoB(k, current);
        r = r * p->lattice_parameter;
        point3D_periodicBoundaryTransform(k, p);
        rdf * r_data = rdf_radius_retrive(data, r);
        int ngbrIndex = point3D_point3DtoIndexFCC(k, p);

        if (a[ngbrIndex] == a[index])
        {
            if (a[ngbrIndex] == 0)
            {
                energy1                 = energy1 + r_data->p11;
                energy2                 = energy2 + r_data->p12;
                chargeDen       = chargeDen + r_data->eDen1;
            }
            else
            {
                energy1                 = energy1 + r_data->p22;
                energy2                 = energy2 + r_data->p12;
                chargeDen       = chargeDen + r_data->eDen2;
            }
        }
        else
        {
            if (a[ngbrIndex] == 0)
            {
                energy1                 = energy1 + r_data->p12;
                energy2                 = energy2 + r_data->p22;
                chargeDen       = chargeDen + r_data->eDen1;
            }
            else
            {
                energy1                 = energy1 + r_data->p12;
                energy2                 = energy2 + r_data->p11;
                chargeDen       = chargeDen + r_data->eDen2;
            }
        }
        free(r_data);
        free(k);
    }

    eDen_df * embeddingEnergy = eDen_df_charge_density_retrive(data, chargeDen);

    if (a[index] == 0)
    {
        energy1 = energy1 + embeddingEnergy->embed1;
        energy2 = energy2 + embeddingEnergy->embed2;
    }
    else
    {
        energy1 = energy1 + embeddingEnergy->embed2;
        energy2 = energy2 + embeddingEnergy->embed1;
    }
    free(embeddingEnergy);
    free(current);
    return energy2 - energy1;
}

/**
 * Function to calculate the energy required to swap an atom at a given index.
 * If say that the atom was initally atom1 then on swapping it would be atom2.
 *
 * @param  index index of ATOM array
 * @param  a     ATOM array
 * @param  data  EAM potential
 * @param  p     parameters of the system
 * @param  ngbrs list of neighbour atoms.
 * @return       change in energy when atom at the given index is swapped
 */
double energyToSwap_fast(int index, int * a, parameter * p, Sn_fcc * ngbrs)
{
    point3D * current = point3D_indexToPoint3D_fcc(index, p);
    int Ni[3] = {0};
    for (size_t n = 0; n < 3; n++)
    {
        for (int i = 0; i < ngbrs->indices[n]; i++)
        {
            point3D * k = point3D_addVectors(current, &(ngbrs->s1n[i]));

            point3D_periodicBoundaryTransform(k, p);
            int ngbrIndex = point3D_point3DtoIndexFCC(k, p);
            if (a[ngbrIndex] == 1)
            {
                Ni[n]++;
            }
            free(k);
        }
    }
    if (a[index] == 0)
    {
        return energyTableInstantLookup[1][Ni[0]][Ni[1]][Ni[2]] - \
                energyTableInstantLookup[0][Ni[0]][Ni[1]][Ni[2]];
    }
    else
    {
        return energyTableInstantLookup[0][Ni[0]][Ni[1]][Ni[2]] - \
                energyTableInstantLookup[1][Ni[0]][Ni[1]][Ni[2]];
    }
}


/**
 * Computes and reports change in energy at every index of the ATOM array
 * @param a            ATOM array
 * @param data         EAM potential data
 * @param p            simulation parameters
 * @param ngbrs        neighbours list
 * @return             pointer to created change in energy matrix
 */
double * deltaEnergyMatrix(int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    double * energyMatrix = (double *) malloc(sizeof(double) * p->no_of_atoms);

    for (int i = 0; i < p->no_of_atoms; i++)
    {
        (energyMatrix)[i] = energyToSwap(i, a, data, p, ngbrs);
    }
    return energyMatrix;
}

/**
 * Creates a lookUpTable for the energy from EAM data so that no extra look up
 * operations are needed for radius fields.
 * @param  data  EAM potential table
 * @param  p     simulation parameters
 * @param  ngbrs FCC nearest neighbours upto 7th
 * @return       pointer to newly created lookup table
 */
lookUpTable * createLookUpTable(binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    int a[] = { 0, 12, 18, 42, 54, 78, 86, };
    lookUpTable * new = malloc(sizeof(lookUpTable));

    for (size_t i = 0; i < 7; i++)
    {
        rdf * temp1;
        double radius = p->lattice_parameter * point3D_magnitude(&(ngbrs->s1n[a[i]]));
        temp1 = rdf_radius_retrive(data, radius);
        new->pairwisePotential_Sn[i] = *(temp1);
        free(temp1);
    }
    return new;
}

void printLookUpTable(lookUpTable * p)
{
    printf("Radius\t\t\tPairwise Potential\t\t\tElectron Denisty\n");
    printf("\t\t  Al-Al \t  Al-Ni \t  Ni-Ni \t  Al  \t\t Ni\n");
    printf("-----------------------------------------------------------------------------------------\n");
    for (size_t i = 0; i < 7; i++)
    {
        printf("%f\t%f\t%f\t%f\t%f\t%f\n", p->pairwisePotential_Sn[i].radius   \
               , p->pairwisePotential_Sn[i].p11      \
               , p->pairwisePotential_Sn[i].p12      \
               , p->pairwisePotential_Sn[i].p22      \
               , p->pairwisePotential_Sn[i].eDen1    \
               , p->pairwisePotential_Sn[i].eDen2);
    }
}
/**
 * Fills the global variable energyTableInstantLookup array with appropriate values
 * @warning This function depends on lookUpTable so it must be not be NULL
 * @param data lookUpTable to read potentials from
 * @param pot  EAM potential data to compute embedding energy.
 * @callgraph
 */
void buildInstantEnergyLookup(lookUpTable *data, binEAMpot *pot)
{
    int i0 = 1;
    for (size_t i1 = 0; i1 < 13; i1++) {
        for (size_t i2 = 0; i2 < 7; i2++) {
            for (size_t i3 = 0; i3 < 25; i3++) {
                energyTableInstantLookup[i0][i1][i2][i3] = i1 * data->pairwisePotential_Sn[0].p22 + \
                                                    (12 - i1) * data->pairwisePotential_Sn[0].p12 + \
                                                           i2 * data->pairwisePotential_Sn[1].p22 + \
                                                     (6 - i2) * data->pairwisePotential_Sn[1].p12 + \
                                                           i3 * data->pairwisePotential_Sn[2].p22 + \
                                                    (24 - i3) * data->pairwisePotential_Sn[2].p12;
                double temp = i1 * data->pairwisePotential_Sn[0].eDen2 + (12 - i1) * data->pairwisePotential_Sn[0].eDen1 + \
                       i2 * data->pairwisePotential_Sn[1].eDen2 + (6 - i2) * data->pairwisePotential_Sn[1].eDen1 + \
                       i3 * data->pairwisePotential_Sn[2].eDen2 + (12 - i3) * data->pairwisePotential_Sn[2].eDen1;
                eDen_df *embeddingE = eDen_df_charge_density_retrive(pot, temp);
                energyTableInstantLookup[i0][i1][i2][i3] += embeddingE->embed2;
            }
        }
    }
    i0 = 0;
    for (size_t i1 = 0; i1 < 13; i1++) {
        for (size_t i2 = 0; i2 < 7; i2++) {
            for (size_t i3 = 0; i3 < 25; i3++) {
                energyTableInstantLookup[i0][i1][i2][i3] = i1 * data->pairwisePotential_Sn[0].p12 + \
                                                    (12 - i1) * data->pairwisePotential_Sn[0].p11 + \
                                                           i2 * data->pairwisePotential_Sn[1].p12 + \
                                                     (6 - i2) * data->pairwisePotential_Sn[1].p11 + \
                                                           i3 * data->pairwisePotential_Sn[2].p12 + \
                                                    (24 - i3) * data->pairwisePotential_Sn[2].p11;
                double temp = i1 * data->pairwisePotential_Sn[0].eDen2 + (12 - i1) * data->pairwisePotential_Sn[0].eDen1 + \
                       i2 * data->pairwisePotential_Sn[1].eDen2 + (6 - i2) * data->pairwisePotential_Sn[1].eDen1 + \
                       i3 * data->pairwisePotential_Sn[2].eDen2 + (12 - i3) * data->pairwisePotential_Sn[2].eDen1;
                eDen_df *embeddingE = eDen_df_charge_density_retrive(pot, temp);
                energyTableInstantLookup[i0][i1][i2][i3] += embeddingE->embed1; }
        }
    }
}
/**
 * Gives concentration of atom1 in a array
 * @param  test input of atomic array
 * @param  p    simulation parameter
 * @return      fractional concentration of atom1
 */
double avgConcentrationAtom1(ATOM *test, parameter *p)
{
    int count = 0;
    for (size_t i = 0; i < p->no_of_atoms; i++) {
        if(test[i] == 0)
            count++;
    }
    return (double)count/(double)p->no_of_atoms;
}

/**
 * Gives number of atom1 in a array
 * @param  test input of atomic array
 * @param  p    simulation parameter
 * @return      fractional concentration of atom1
 */

int atomsType1(ATOM *test, parameter *p)
{
    int count = 0;
    for (size_t i = 0; i < p->no_of_atoms; i++) {
        if(test[i] == 0)
            count++;
    }
    return count;
}
