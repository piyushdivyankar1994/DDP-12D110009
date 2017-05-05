/*!
   \file parameters.c
   \brief Source file for functions in parameters.h
   \author Piyush Divyankar
   \date 01/09/2016
*/

#include "parameters.h"

/**
 * Prints the ::parameter data to STDOUT.
 * @param a pointer to parameter
 */
void print_parameters(parameter * a)
{
    printf("# Size in x-direction           = %d\n", a->Nx);
    printf("# Size in y-direction           = %d\n", a->Ny);
    printf("# Size in z-direction           = %d\n", a->Nz);
    printf("# Lattice Parameter             = %le\n", a->lattice_parameter);
    printf("# No. of MonteCarlo Simulations = %d\n", a->N_MCS);
    printf("# Temperature                   = %le\n", a->temperature);
    printf("# Atoms per site                = %d\n", a->atoms_per_site);
    printf("# Totol No of atoms             = %d\n", a->no_of_atoms);
}

/**
 * Reads a parameter from given file name.
 * @param  filename string that contains the file name.
 * @return          pointer to newly created parameter struct
 */
parameter * new_parameters(char * fileName)
{
    parameter * new = (parameter *) malloc(sizeof(parameter));
    FILE * fp;

    fp = fopen(fileName, "r");
    int N;

    fscanf(fp, "%d", &N);
    new->Nx = N;
    printf("%d\n", N);

    fscanf(fp, "%d", &(new->Ny));
    fscanf(fp, "%d", &(new->Nz));
    fscanf(fp, "%le", &(new->lattice_parameter));
    fscanf(fp, "%d", &(new->N_MCS));
    fscanf(fp, "%le", &(new->temperature));
    fscanf(fp, "%d", &(new->atoms_per_site));

    new->nearestNeighbours[0] = 12;
    new->nearestNeighbours[1] = 6;
    new->nearestNeighbours[2] = 24;
    new->nearestNeighbours[3] = 12;
    new->nearestNeighbours[4] = 24;
    new->nearestNeighbours[5] = 8;
    new->nearestNeighbours[6] = 48;
    new->no_of_atoms = new->Nx * new->Ny * new->Nz * new->atoms_per_site;
    strcpy(new->fileName, fileName);
    fclose(fp);
    return new;
}

/**
 * Creates a ::parameter with default values as follows truct as follows.
 * Size in x-direction           = ::parameter.Nx = 10 <BR>
 * Size in y-direction           = ::parameter.Ny = 10 <BR>
 * Size in z-direction           = ::parameter.Nz = 10 <BR>
 * Lattice Parameter             = ::parameter.lattice_parameter = 4.0 A <BR>
 * No. of MonteCarlo Simulations = ::parameter.N_MCS             = 1000  <BR>
 * Temperature                   = ::parameter.temperature       = 1000  <BR>
 * Atoms per site                = ::parameter.atoms_per_site    = 4     <BR>
 * Totol No of atoms             = ::parameter.no_of_atoms       = 4000  <BR>
 * @return pointer to newly created parameter struct1
 */
parameter * _defaultFCCparameter()
{
    parameter * new = (parameter *) malloc(sizeof(parameter));

    new->N_MCS = 100;
    new->Nx    = 20;
    new->Ny    = 20;
    new->Nz    = 20;
    new->atoms_per_site    = 4;
    new->lattice_parameter = 4.00;
    new->temperature       = 1000;
    new->nearestNeighbours[0] = 12;
    new->nearestNeighbours[1] = 6;
    new->nearestNeighbours[2] = 24;
    new->nearestNeighbours[3] = 12;
    new->nearestNeighbours[4] = 24;
    new->nearestNeighbours[5] = 8;
    new->nearestNeighbours[6] = 48;
    new->no_of_atoms          = new->Nx * new->Ny * new->Nz * new->atoms_per_site;
    strcpy(new->fileName, "defaultFCC.param");
    return new;
}

parameter * _defaultBCCparameter()
{
    parameter * new = (parameter *) malloc(sizeof(parameter));

    new->N_MCS = 100;
    new->Nx    = 10;
    new->Ny    = 10;
    new->Nz    = 10;
    new->atoms_per_site    = 2;
    new->lattice_parameter = 4.00;
    new->temperature       = 1000;
    new->nearestNeighbours[0] = 8;
    new->nearestNeighbours[1] = 6;
    new->nearestNeighbours[2] = 24;
    new->nearestNeighbours[3] = 12;
    new->nearestNeighbours[4] = 24;
    new->nearestNeighbours[5] = 8;
    new->nearestNeighbours[6] = 48;
    new->no_of_atoms          = new->Nx * new->Ny * new->Nz * new->atoms_per_site;
    strcpy(new->fileName, "defaultBCC.param");
    return new;
}


/**
 * Creates a ::parameter structure from STDIN and saves it to a .parameter FILE
 * @return pointer to new parameter obejct.
 */
parameter * createParameterFileFromInput()
{
    parameter * new = (parameter *) malloc(sizeof(parameter));

    new = _defaultFCCparameter();

    printf("\nEnter # of MonteCarlo simulations per atom = ");
    scanf("%d", &(new->N_MCS));

    printf("\nEnter Lattice size in x-direction = ");
    scanf("%d", &(new->Nx));

    printf("\nEnter Lattice size in y-direction = ");
    scanf("%d", &(new->Ny));

    printf("\nEnter Lattice size in z-direction = ");
    scanf("%d", &(new->Nz));

    printf("\nEnter atoms per site = ");
    scanf("%d", &(new->atoms_per_site));

    printf("\nEnter lattice parameter = ");
    scanf("%le", &(new->lattice_parameter));

    printf("\nEnter simulation Temperature");
    scanf("%le", &(new->temperature));

    new->no_of_atoms = new->Nx * new->Ny * new->Nz * new->atoms_per_site;

    new->nearestNeighbours[0] = 12;
    new->nearestNeighbours[1] = 6;
    new->nearestNeighbours[2] = 24;
    new->nearestNeighbours[3] = 12;
    new->nearestNeighbours[4] = 24;
    new->nearestNeighbours[5] = 8;
    new->nearestNeighbours[6] = 48;

    /** Name of the FILE must be less than 50 characters */
    char fileName[50];

    printf("\nEnter a filename for list of parameters");
    scanf("%s", fileName);
    strcat(fileName, ".param");
    strcpy(new->fileName, fileName);
    FILE * fp = fopen(fileName, "w");
    fwrite(new, sizeof(parameter), 1, fp);

    fclose(fp);

    return new;
}

/**
 * Reads input from a file that is created by ::createParameterFileFromInput(...)
 * @param fileName name of a valid file.
 * @return pointer to newly created parameter
 * @warning Code will return NULL for bad fileName and files whose size doesn't
 * match that of parameter.
 */
parameter * parameterReadFromFile(char fileName[50])
{
    FILE * fp = fopen(fileName, "r");

    // / FIXME:0 Using standard C library throw an error here.
    if (fp == NULL)
    {
        return NULL;
    }

    parameter * new = _defaultFCCparameter();

    fseek(fp, 0L, SEEK_END);
    if (sizeof(parameter) != ftell(fp))
    {
        free(new);
        fclose(fp);
        return NULL;
    }

    fread(new, sizeof(parameter), 1, fp);
    fclose(fp);
    return new;
}

/**
 * Creates default ::parameter file named defaultFCC.param in working directory.
 */
void parameterDefaultFile()
{
    parameter * new = _defaultFCCparameter();

    printf("\nCreating a default file with\n");
    print_parameters(new);

    FILE * fp = fopen(new->fileName, "w");
    fwrite(new, sizeof(parameter), 1, fp);
    fclose(fp);
    free(new);
    return;
}
/**
 * Takes a ::parameter object and writes it to a file.
 * @param out      pointer to parameter
 */
void parameterWriteToFile(parameter * out)
{
    FILE * fp = fopen(out->fileName, "w");

    fwrite(out, sizeof(parameter), 1, fp);
    fclose(fp);
}
