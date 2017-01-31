#include "analysis.h"

/**
 * Calculates total bond energy in a atomic matrix.
 * @param  a     ::ATOM array(poistion of various atoms in lattice)
 * @param  data  ::binaryEAMpotential for computing bond energies.
 * @param  p     ::parameter for the given lattice.
 * @param  ngbrs ::neighbours_fcc containing neighbours relative to (0, 0, 0) in fcc crystal
 * @return       net bond energy
 * @pre          It is assumed that atomic matrix, binaryEAMpotential, parameter,
 *               and neighbours_fcc are consistent. No checks are made about their
 *               consistency.
 */
double analysis_totalEnergy(int * a, binEAMpot * data, parameter * p, Sn_fcc * ngbrs)
{
    double ret_val = 0;

    for (size_t i = 0; i < p->no_of_atoms; i++)
    {
        ret_val = ret_val + energyAtIndexFCC_ver2(i, a, data, p, ngbrs, 12);
    }

    return ret_val;
}

/**
 * Counts total number of ordered phase sites in atomic matrix
 * @param  a     ::ATOM array(poistion of various atoms in lattice)
 * @param  p     ::parameter for the given lattice.
 * @param  ngbrs ::neighbours_fcc containing neighbours relative to (0, 0, 0) in fcc crystal
 * @return       integer number of ordered phase counts.
 */
int orderedPhaseCount(int * a, parameter * p, Sn_fcc * ngbrs)
{
    printf("oregnioegn\n");
    int count = 0;
    int j;
    for (size_t i = 0; i < p->no_of_atoms; i++)
    {
        point3D * new = point3D_indexToPoint3D_fcc(i, p);
        for (j = 0; j < 12; j++)
        {
            point3D * test = point3D_addVectors(new, &(ngbrs->s1n[j]));
            int t = point3D_point3DtoIndexFCC(test, p);
            if (a[t] == a[i])
            {
                break;
            }
            free(test);
        }
        if (j == 12)
        {
            count++;
        }
        free(new);
    }
    printf("oregnioegn\n");
    return count;
}

/**
 * Counts total number of anti-ordered phase sites in atomic matrix
 * @param  a     ::ATOM array(poistion of various atoms in lattice)
 * @param  p     ::parameter for the given lattice.
 * @param  ngbrs ::neighbours_fcc containing neighbours relative to (0, 0, 0) in fcc crystal
 * @return       integer number of anti-ordered phase counts.
 */
int antiOrderedPhaseCount(int * a, parameter * p, Sn_fcc * fcc)
{
    int count = 0;
    int j;
    int checkSite = 0;

    for (size_t i = 0; i < p->no_of_atoms; i++)
    {
        point3D * new = point3D_indexToPoint3D_fcc(i, p);
        checkSite = 0;
        for (j = 0; j < 12; j++)
        {
            point3D * test = point3D_addVectors(new, &(fcc->s1n[j]));
            int t = point3D_point3DtoIndexFCC(test, p);
            if (a[t] != a[i])
            {
                checkSite++;
            }
        }
        if (checkSite == 11)
        {
            printf("flag %d \n", j);
            count++;
        }
    }
    return count;
}

/**
 * Generates a random atomic matrix.
 * @param p               ::parameter of matrix.
 * @param ouput_file_name output file name(written to /"output_file_name".crystal.fcc)
 * @param seed            [description]
 */
void randomMatrixGeneratorFCC(parameter * p, char * ouput_file_name, unsigned long int seed, double concentration)
{
    const gsl_rng_type * T;
    gsl_rng * r;

    //strcat(ouput_file_name, ".crystal.fcc");
    FILE * fp = fopen(ouput_file_name, "w");
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    for (size_t i = 0; i < p->no_of_atoms; i++)
    {
        double u = gsl_rng_uniform(r);
        int a;
        if (u > concentration)
        {
            a = 1;
        }
        else
        {
            a = 0;
        }
        if (fwrite(&a, sizeof(int), 1, fp) == 1)
        {
            continue;
        }
        else
        {
            i--;
        }
    }
    char str[1000];
    strcpy(str, ouput_file_name);
    strcat(str, ".random.gen.param");
    FILE * fr = fopen(str, "w");
    gsl_rng_fwrite(fr, r);
    fclose(fp);
    gsl_rng_free(r);
    fclose(fr);
}

/**
 * Counts total number of atoms in a file
 * @param  fp file of the form *.crystal.fcc
 * @return    number of atoms
 */
int totalAtomsInFile(FILE * fp)
{
    fseek(fp, 0L, SEEK_END);
    int no_of_atoms = ftell(fp) / sizeof(int);
    fseek(fp, 0L, SEEK_SET);
    return no_of_atoms;
}

/**
 * Calculates total energy using energy look up table
 * @param  a     Atomic matrix
 * @param  p     Simulation parameters
 * @param  ngbrs List of neighbours
 * @return       total bond energy of the lattice
 */
double totalEnergyQuick(ATOM * a, parameter * p, Sn_fcc * ngbrs)
{
    double ret_val = 0;
    for (size_t i = 0; i < p->no_of_atoms; i++) {
        ret_val += energyAtIndexFCC_fast(i, a, p, ngbrs);
    }
    return ret_val;
}
void L12lattice(parameter *p, char *filename, unsigned long int seed)
{
    //strcat(ouput_file_name, ".crystal.fcc");
    FILE * fp = fopen(filename, "w");
    for (size_t i = 0; i < p->no_of_atoms; i++)
    {
        int a;
        if (i % 4 != 0)
        {
            a = 1;
        }
        else
        {
            a = 0;
        }
        if (fwrite(&a, sizeof(int), 1, fp) == 1)
        {
            continue;
        }
        else
        {
            i--;
        }
    }

    fclose(fp);
    gsl_rng_free(r);
}
