#include "test.h"
#include "simulation.h"
// / FUTURE_CHANGES:40 Think of a simulationData type of struct.
// / FUTURE_CHANGES:50 Two phase equlibria simulation.
// / FIXME:10 Analysis code for two phase equilibria.
// / FUTURE_CHANGES:60 chemicalPotentialAtIndex(...) is incorrectly evaluated correct it
// / TODO:10 Place error messages at places where file names need to be passed among functions.
// / TODO:30 documentation
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
    test_parametersInputOutput();
    return 0;
}
