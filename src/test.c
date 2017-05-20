/*!
   \file test.c
   \brief Source file for test functions
   \author Piyush Divyankar
   \date 01/09/2016
 */
#include "test.h"

void test_new_parameters()
{
    parameter * new = NULL;

    new = new_parameters("parameter1.txt");
    print_parameters(new);
    free(new);
}

void test_dataRetrival()
{
    double r = 4.123;
    binEAMpot * data = NULL;
    rdf * nn = (rdf *) malloc(sizeof(rdf));
    eDen_df * a = (eDen_df *) malloc(sizeof(eDen_df));

    data = eam_data_read("file_list.txt", "Al", "Ni");
    nn = rdf_radius_retrive(data, r);
    printf("%d\n", nn->index);
    printf("%le\n", nn->p11);
    printf("%le\n", nn->p12);
    printf("%le\n", nn->p22);
    printf("%le\n", nn->eDen1);
    printf("%le\n", nn->eDen2);

    a = eDen_df_charge_density_retrive(data, 0.01);
    printf("\n");
    printf("%le\n", a->eDen);
    printf("%le\n", a->embed1);
    printf("%le\n", a->embed2);
    free(data);
}

void test_eam_data_read()
{
    binEAMpot * data = NULL;

    data = eam_data_read("/home/piyushdivyankar/Desktop/DDP-12D110009/EAM_Ni_Al/file_list.txt", "Al", "Ni");

    printf("%d\n", data->no_of_files);
    printf("%s\n", data->atom1);
    printf("%s\n", data->atom2);
    free(data);
}

void test_point3D_indexToPoint3D_fcc()
{
    parameter * p = _defaultFCCparameter();

    printf("Generating 10 random numbers from 0 to 3999 and displaying corresponding points\n");
    for (int i = 0; i < 15; i++)
    {
        int index = rand() % p->no_of_atoms;
        printf("%d\t=\t", index);
        point3D_dispPoint(point3D_indexToPoint3D_fcc(index, p));
    }
}

void test_point3D()
{
    point3D * a1 = point3D_origin();

    point3D_dispPoint(a1);
    a1 = point3D_newPoint(1, 1.5, 4.3);
    point3D_dispPoint(a1);
    point3D * a2 = point3D_newPoint(1.1, 3.4, 5.3);
    point3D_dispPoint(a2);
    point3D * sum = point3D_addVectors(a1, a2);
    point3D_dispPoint(sum);
    point3D * diff = point3D_subtractVectors(a1, a2);
    point3D_dispPoint(diff);
    printf("Distance Test from point3D_origin\n");
    a1 = point3D_newPoint(1, 1, 1);
    printf("%f\n", point3D_magnitude(a1));
    a1 = point3D_newPoint(1, 1.5, 1);
    printf("%f\n", point3D_magnitude(a1));
    a1 = point3D_newPoint(1.4, 1, 1);
    printf("%f\n", point3D_magnitude(a1));
    a1 = point3D_newPoint(0, 0, 5);
    printf("%f\n", point3D_magnitude(a1));

    printf("Distance between two points\n");
    for (int i = 0; i < 3; i++)
    {
        a1 = point3D_newPoint(rand() / (float) RAND_MAX, rand() / (float) RAND_MAX, rand() / (float) RAND_MAX);
        a2 = point3D_newPoint(rand() / (float) RAND_MAX, rand() / (float) RAND_MAX, rand() / (float) RAND_MAX);
        point3D_dispPoint(a1);
        point3D_dispPoint(a2);
        printf("Distance = %f\n\n", point3D_distAtoB(a1, a2));

    }

    return;
}

void test_defaultFCCparameter()
{
    parameter * new = _defaultFCCparameter();

    print_parameters(new);
}


void test_Sn_fcc_readNeighbours_fromFile()
{
    Sn_fcc * new = Sn_fcc_readNeighbours_fromFile("file_list_neighbours.txt");

    // print_Neighbours(new);
    printf("Everything is continuoulsy stored so if i goes from 0 \
	to sum(indicies) then s1n to s7n can be accessed as s1n[i] \
	this is demonstrated by printing 13 points\n "    );
    for (int i = 0; i < 13; i++)
    {
        point3D_dispPoint(&(new->s1n[i]));
    }
}


void test_AtomicMatrixRead()
{
    parameter * p = _defaultFCCparameter();
    int * a = atomicMatrixRead("out.txt", p);

    print_AtomicMatrix(a, 20, 50);
}


void test_energyAtIndexFCC()
{
    int index = 100;
    binEAMpot * data = NULL;

    data = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    // test_AtomicMatrixRead();
    Sn_fcc * new = _defaultFCCNeighbours();
    double e = energyAtIndexFCC(index, a, data, p, new);
    printf("at index = %d energy = %f\n", index, e);
    free(data);
    free(p);
    free(a);
    free(new);
}

void test_point3D_point3DtoIndexFCCFCC()
{
    parameter * p = _defaultFCCparameter();
    point3D * p1 = point3D_newPoint(2, 0, 0);
    point3D * p2 = point3D_newPoint(0.5, 0.5, 1);
    point3D * p3 = point3D_newPoint(2.0, 0.5, 0.5);
    point3D * p4 = point3D_newPoint(0.5, 1, 0.5);

    printf("%d\t", point3D_point3DtoIndexFCC(p1, p));
    point3D_dispPoint(p1);
    printf("%d\t", point3D_point3DtoIndexFCC(p2, p));
    point3D_dispPoint(p2);
    printf("%d\t", point3D_point3DtoIndexFCC(p3, p));
    point3D_dispPoint(p3);
    printf("%d\t", point3D_point3DtoIndexFCC(p4, p));
    point3D_dispPoint(p4);
}

void test_point3D_periodicBoundaryTransform()
{
    point3D * new = point3D_newPoint(-10.5, 11.5, 16.5);
    parameter * p = _defaultFCCparameter();

    printf("test function for periodic boundary transform\n");
    point3D_dispPoint(new);
    // point3D* test =
    point3D_periodicBoundaryTransform(new, p);
    point3D_dispPoint(new);
}

// / TODO find memory leak in this
void test_energyInMatrix()
{
    binEAMpot * data = NULL;

    data = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    // test_AtomicMatrixRead();
    Sn_fcc * new = _defaultFCCNeighbours();

    // / Example:
    double * energyMap = energyInMatrix(a, data, p, new);
    printEnergyMap(energyMap, 100, 200);
    free(p);
    free(a);
    free(new);
    free(energyMap);

}

// / DONE:10 Memory leaks in this function.
// / Occured due to ::energyToSwap function now fixed
void test_deltaEnergyMatrix()
{
    binEAMpot * data = NULL;

    data = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    // test_AtomicMatrixRead();
    Sn_fcc * new = _defaultFCCNeighbours();

    // / Example:
    double * energyMap = deltaEnergyMatrix(a, data, p, new);
    printEnergyMap(energyMap, 100, 200);
    free(data);
    free(p);
    free(a);
    free(new);
    free(energyMap);
}

/*!
   \brief Tests chemicalPotentialAtIndex function. Takes in a matrix and calculates
   chemical potential at 100-199 indicies.
 */

void test_chemicalPotentialAtIndex()
{
    binEAMpot * data = NULL;

    data = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    Sn_fcc * new = _defaultFCCNeighbours();

    for (int i = 100; i < 199; i++)
    {
        printf("%le\n", chemicalPotentialAtIndex(i, a, data, p, new));
    }
}

void test_analysis_totalEnergy()
{
    binEAMpot * data = NULL;

    data = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    Sn_fcc * new = _defaultFCCNeighbours();

    printf("%le\n", analysis_totalEnergy(a, data, p, new));
}

void test_orderedPhaseCount()
{
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    Sn_fcc * fcc = _defaultFCCNeighbours();
    int count = orderedPhaseCount(a, p, fcc);

    printf("ekneiogne\n");
    printf("Order Phase sites in given lattice = %d\n", count);
    free(a);
    free(fcc);
    free(p);
    return;
}

void test_antiOrderedPhaseCount()
{
    parameter * p = _defaultFCCparameter();
    Sn_fcc * fcc = _defaultFCCNeighbours();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    int count = antiOrderedPhaseCount(a, p, fcc);

    printf("Anti-Order Phase sites in given lattice = %d\n", count);
}

void test_randomMatrixGeneratorFCC()
{
    parameter * p = _defaultFCCparameter();
    char fileName[] = "input";

    randomMatrixGeneratorFCC(p, fileName, 12434533, 0.5);
}

void test_readCrystalFileFCC()
{
    char fileName[] = "input.crystal.fcc";
    int * a_file = readCrystalFileFCC(fileName);

    print_AtomicMatrix(a_file, 0, 15);
}

void test_point3D_origin()
{
    point3D * new = point3D_origin();

    point3D_dispPoint(new);
    free(new);
}

void test_point3D_newPoint()
{
    point3D * new = point3D_newPoint(1.5, 1.2, 1.4);

    point3D_dispPoint(new);
    free(new);
}

// / point3D_addVectors(...) should be made equal to one of the argument.
// / It has been shown here in line 4. Instruction on line 4 will cause the
// / memory that new1 was pointing to be unrefrenced as cause a leak.
// / Correct way to do this is to declare a third point3D* and use it to point
// / to the output of the point3D_addVectors function.
void test_point3D_addVectors()
{
    // / Incorrect usage
    point3D * new1 = point3D_newPoint(1, 1, 1);
    point3D * new2 = point3D_newPoint(1, 3, 2);

    new1 = point3D_addVectors(new1, new2);
    free(new1);
    free(new2);

    // / Correct usage.
    new1 = point3D_newPoint(1, 1, 1);
    new2 = point3D_newPoint(1, 3, 2);
    point3D * new3 = point3D_addVectors(new1, new2);
    free(new1);
    free(new2);
    free(new3);
}

/**
 * Tests the functions ::createParameterFileFromInput, ::parameterReadFromFile,
 * ::parameterWriteToFile, and ::parameterDefaultFile
 */
void test_parametersInputOutput()
{
    printf("\nThis is TEST function test_parametersInputOutput() refer documentation for details\n");
    printf("Begining TEST\n");
    printf("This is a STDIN TEST please follow on screen Instructions\n");
    parameter * new = _defaultFCCparameter();
    printf("All ok\n");
    printf("Writing to file\n");
    parameterWriteToFile(new);
    parameter * check = NULL;
    // /DONE:30 Segmentation fault
    check = parameterReadFromFile(new->fileName);
    if (check == NULL)
    {
        printf("flag\n");
        printf("%s\n", new->fileName);
        return;
    }
    print_parameters(check);
    print_parameters(new);
    printf("If two exactly same streams of data appear then its fine\n");
    remove(new->fileName);
    free(new);
    free(check);
    parameterDefaultFile();
    new = parameterReadFromFile(PARAM_FILE_DEF_NAME);
    free(new);
}

void test_createLookUpTable()
{
    binEAMpot * data     = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * defFCC   = _defaultFCCNeighbours();
    parameter * simParam = parameterReadFromFile("parametersSim1.param");

    simParam->lattice_parameter = 3.56;
    lookUpTable * table    = createLookUpTable(data, simParam, defFCC);
    printLookUpTable(table);
}

/**
 * Following function generates data and has been used to show that taking upto
 * 4th nearest neighbours is enough accuracy for calculation.
 */
void test_nearest_neighbours()
{
    binEAMpot * data     = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * defFCC   = _defaultFCCNeighbours();
    parameter * simParam = parameterReadFromFile("parametersSim1.param");

    simParam->lattice_parameter = 3.56;


    randomMatrixGeneratorFCC(simParam, "input", 1232342, 0.5);
    ATOM * test_input = readCrystalFileFCC("input");
    int a[] = { 12, 18, 42, 54, 78, 86, 134 };
    double * S7energy = energyInMatrix(test_input, data, simParam, defFCC);
    for (size_t j = 0; j < 7; j++)
    {
        double * S4energy = energyInMatrix_ver2(test_input, data, simParam, defFCC, a[j]);

        /** Taking statistics */
        double totalE7 = 0;
        double totalE4 = 0;
        double absoulte_mean = 0;
        double variance      = 0;
        for (size_t i = 0; i < simParam->no_of_atoms; i++)
        {
            absoulte_mean += (fabs(S7energy[i] - S4energy[i]) / fabs(S7energy[i]));
            totalE4 += S4energy[i];
            totalE7 += S7energy[i];
        }
        absoulte_mean /= 32000;
        for (size_t i = 0; i < simParam->no_of_atoms; i++)
        {
            variance += pow((absoulte_mean - (fabs(S7energy[i] - S4energy[i]) / fabs(S7energy[i]))), 2);
        }
        variance /= 32000;
        printf("%le, %le\n", absoulte_mean, sqrt(variance));
    }
    // printf("Total E4 = %le Total E7 = %le\n", totalE4, totalE7);
    free(data);
    free(defFCC);
    free(simParam);
}

void test_buildInstantEnergyLookup()
{
    binEAMpot * data     = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * defFCC   = _defaultFCCNeighbours();
    parameter * simParam = parameterReadFromFile("parametersSim1.param");

    simParam->lattice_parameter = 3.56;

    lookUpTable * t = createLookUpTable(data, simParam, defFCC);
    buildInstantEnergyLookup(t, data);
    for (size_t i0 = 0; i0 < 2; i0++)
    {
        for (size_t i1 = 0; i1 < 13; i1++)
        {
            for (size_t i2 = 0; i2 < 7; i2++)
            {
                for (size_t i3 = 0; i3 < 25; i3++)
                {
                    printf("%le\n", energyTableInstantLookup[i0][i1][i2][i3]);
                }
            }
        }
    }
}

void test_energyAtIndexFCC_fast()
{
    binEAMpot * data     = eam_data_read("./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * defFCC   = _defaultFCCNeighbours();
    parameter * simParam = parameterReadFromFile("parametersSim1.param");

    simParam->lattice_parameter = 3.56;

    lookUpTable * t = createLookUpTable(data, simParam, defFCC);
    buildInstantEnergyLookup(t, data);

    randomMatrixGeneratorFCC(simParam, "input", 1232342, 0.5);
    ATOM * test_input = readCrystalFileFCC("input");

    for (size_t i = 0; i < 32000; i++)
    {
        printf("%le\n", energyAtIndexFCC_fast(i, test_input, simParam, defFCC));
    }
}

void test_point3D_neighbourIndexTableFCC()
{
    Sn_fcc * defFCC   = _defaultFCCNeighbours();
    parameter * simParam = parameterReadFromFile("parametersSim1.param");
    int * ngbrTable = point3D_neighbourIndexTable_FCC(simParam, defFCC, 3);

    for (int i = 0; i < simParam->no_of_atoms; i++)
    {
        for (int j = 0; j < 42; j++)
        {
            printf("%d, ", ngbrTable[i * 42 + j]);
        }
        printf("\n");
    }
}

void test_concentrationFunctions()
{
    ATOM * inputMatrix = readCrystalFileFCC("inputCrystalFiles/input");
    Sn_fcc * defFCC   = _defaultFCCNeighbours();
    parameter * simParam = parameterReadFromFile("parametersSim1.param");
    int * ngbrTable = point3D_neighbourIndexTable_FCC(simParam, defFCC, 3);
    for(int i = 0; i < simParam->no_of_atoms; i += 100) {
        printf("%0.3f\n", concentrationAtIndex(inputMatrix, ngbrTable , i, 42));
       }
    double * table = createConcentrationTable(inputMatrix, simParam, ngbrTable, 42);

    // printConcentrationTable(table, 0, 10);
    for (int i = 0; i < 42; i++)
    {
        printf("%0.5f\n", table[ngbrTable[i]]);
    }
    printf("\n");
    updateConcentrationTable(table, ngbrTable, 0, 42, 1);
    for (int i = 0; i < 42; i++)
    {
        printf("%0.5f\n", table[ngbrTable[i]]);
    }
    printf("\n");
}

void test_displayInstantEnergyLookUpTable()
{
    parameter * AlNi_fcc = parameterReadFromFile("parametersSim1.param");
    binEAMpot * potential = NULL;

    potential = eam_data_read("EAM_Ni_Al/file_list.txt", "Al", "Ni");
    Sn_fcc * fccNeighbours = _defaultFCCNeighbours();

    lookUpTable * t = createLookUpTable(potential, AlNi_fcc, fccNeighbours);
    buildInstantEnergyLookup(t, potential);
    for (size_t i0 = 0; i0 < 2; i0++)
    {
        for (size_t i1 = 0; i1 < 13; i1++)
        {
            for (size_t i2 = 0; i2 < 7; i2++)
            {
                for (size_t i3 = 0; i3 < 25; i3++)
                {
                    printf("%le\n", energyTableInstantLookup[i0][i1][i2][i3]);
                }
            }
        }
    }
    return;
}

void test_point3D_indexToPoint3D_bcc() {
  parameter * pBCC = _defaultFCCparameter();
  pBCC->atoms_per_site = 2;
  pBCC->no_of_atoms = 2 * pBCC->Nx * pBCC->Ny * pBCC->Nz;

  printf("Generating 10 random numbers from 0 to 3999 and displaying corresponding points\n");
  for (int i = 0; i < 19; i++)
  {
      int index = rand() % pBCC->no_of_atoms;
      point3D * a = point3D_indexToPoint3D_bcc(index, pBCC);
      int k = point3D_point3DtoIndexBCC(a, pBCC);
      printf("%d\t=\t, [%d] ", index, k);
      point3D_dispPoint(a);
      free(a);
      //point3D_dispPoint(point3D_indexToPoint3D_bcc(i, pBCC));
  }
}

void test_negihbourReading_transformations() {
    parameter * pBCC = _defaultFCCparameter();
    pBCC->atoms_per_site = 2;
    pBCC->no_of_atoms = 2 * pBCC->Nx * pBCC->Ny * pBCC->Nz;
    print_parameters(pBCC);
    Sn_bcc * ngbrUC = readBCCfromFile( "/home/piyush/Desktop/DDP-12D110009/neighbours/bcc/bccNeighbours.txt");
    point3D * a = point3D_origin();
    /*for(int i = 0; i < 14; i++)
    {
        point3D * na = point3D_addVectors(a, &ngbrUC->s1n[i]);
        //point3D_dispPoint(na);
        //point3D_dispPoint(&(ngbrUC->s1n[i]));
        point3D_periodicBoundaryTransform(na, pBCC);
        int k = point3D_point3DtoIndexBCC(na, pBCC);
        printf("## [%d]\t=\t\n", k);
        point3D_dispPoint(na);
        free(na);
    }*/

    int * ngbrTable = point3D_neighbourIndexTable_BCC(pBCC,ngbrUC);
    for(int i = 0; i < 100; i++) {
        printf("[%d]\t-> ", i);
        for(int j = 0; j < 14; j++) {
            printf("[%d], ", ngbrTable[i*14+j]);
        }
        printf("\n");
    }
    free(a);
    free(ngbrTable);
    free(pBCC);
    free(ngbrUC);

}
