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

    eam_data_read(&data, "file_list.txt", "Al", "Ni");
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
}

void test_eam_data_read()
{
    binEAMpot * data = NULL;

    eam_data_read(&data, "file_list.txt", "Al", "Ni");

    printf("%d\n", data->no_of_files);
    printf("%s\n", data->atom1);
    printf("%s\n", data->atom2);
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

    eam_data_read(&data, "./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    // test_AtomicMatrixRead();
    Sn_fcc * new = _defaultFCCNeighbours();
    double e = energyAtIndexFCC(index, a, data, p, new);
    printf("at index = %d energy = %f\n", index, e);

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

void test_energyInMatrix()
{
    binEAMpot * data = NULL;

    eam_data_read(&data, "./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    // test_AtomicMatrixRead();
    Sn_fcc * new = _defaultFCCNeighbours();

    // / Example:
    double * energyMap = NULL;
    energyInMatrix(&energyMap, a, data, p, new);
    printEnergyMap(energyMap, 100, 200);
    free(energyMap);

}

void test_deltaEnergyMatrix()
{
    binEAMpot * data = NULL;

    eam_data_read(&data, "./EAM_Ni_Al/file_list.txt", "Al", "Ni");
    parameter * p = _defaultFCCparameter();
    int * a = readCrystalFileFCC("inputCrystalFiles/input.crystal.fcc");
    // test_AtomicMatrixRead();
    Sn_fcc * new = _defaultFCCNeighbours();

    // / Example:
    double * energyMap = NULL;
    deltaEnergyMatrix(&energyMap, a, data, p, new);
    printEnergyMap(energyMap, 100, 200);
    free(energyMap);
}

/*!
   \brief Tests chemicalPotentialAtIndex function. Takes in a matrix and calculates
   chemical potential at 100-199 indicies.
 */

void test_chemicalPotentialAtIndex()
{
    binEAMpot * data = NULL;

    eam_data_read(&data, "./EAM_Ni_Al/file_list.txt", "Al", "Ni");
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

    eam_data_read(&data, "./EAM_Ni_Al/file_list.txt", "Al", "Ni");
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

    randomMatrixGeneratorFCC(p, fileName, 12434533);
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
    parameter *new = _defaultFCCparameter();
    printf("All ok\n");
    printf("Writing to file\n");
    parameterWriteToFile(new);
    parameter *check = NULL;
    ///DOING:10 Segmentation fault
    check = parameterReadFromFile(new->fileName);
    if (check == NULL) {
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
