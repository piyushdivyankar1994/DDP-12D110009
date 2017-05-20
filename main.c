#include "test.h"
#include "simulation.h"


void thermalExpansion() {
    double temp = 20;
    while( temp < 501) {
        latticeParameterSimulation(rand(), temp);
        temp += 20;
    }
}




int main(int argc, char const * argv[])
{
    // /test_point3D();
    // test_point3D_origin();
    // test_point3D_addVectors();

    // test_energyAtIndexFCC();
     // test_energyInMatrix();
    // test_deltaEnergyMatrix();
    // / FUTURE_CHANGES:0 As of 22 August this fucntion is incorrectly written
    // test_chemicalPotentialAtIndex();

    // test_analysis_totalEnergy();
    // test_orderedPhaseCount();
    /*test_antiOrderedPhaseCount();
       test_randomMatrixGeneratorFCC();
       test_readCrystalFileFCC();
       twoPhaseEquilibriaSimulation(4231332);*/
    //test_parametersInputOutput();
    //test_eam_data_read();
    ////twoPhaseEquilibriaSimulation(123314);
    /*FILE *fp = NULL;
    parameter * new = createParameterFileFromInput();
    fp = fopen(new->fileName, "w");
    fwrite(new, sizeof(parameter), 1, fp);
    fclose(fp);
    free(new);*/
    //est_createLookUpTable();
    //test_buildInstantEnergyLookup();
    //test_energyAtIndexFCC_fast();
    long long int seed_value = 102234;
    //cannonicalEnsemble(seed_value);

    // latticeParameterSimulation(seed_value, 100);
    // test_point3D_neighbourIndexTable();
    //
    // test_concentrationFunctions();
    // semiGrandCanonical_concentration_study(1233);
    // test_displayInstantEnergyLookUpTable();
    //pairwiseConstants();
    // semiGrandCanonical_concentration_study_ljp(1313);
    //ljp_cannonical_order_disorder_transformations(48165532194);


    //test_point3D_periodicBoundaryTransform();
    //test_point3D_indexToPoint3D_bcc();
    //test_negihbourReading_transformations();
    //bccCannonicalBenchmark();
    //bcc_SGCannonicalBenchmark();
    thermalExpansion();
    return 0;
}
