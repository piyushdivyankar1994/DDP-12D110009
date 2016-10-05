#include "test.h"
#include "simulation.h"
// / DOING:0 Two phase equlibria simulation.
// / FUTURE_CHANGES:30 Think of a simulationData type of struct.
// / FIXME:30 Analysis code for two phase equilibria.
// / FUTURE_CHANGES:40 chemicalPotentialAtIndex(...) is incorrectly evaluated correct it
// / TODO:10 Place error messages at places where file names need to be passed among functions.
// / DONE:0 documentation
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
    long long int seed_value = 0;
    cannonicalEnsemble(seed_value);

    //latticeParameterSimulation(32123414);
    return 0;
}
