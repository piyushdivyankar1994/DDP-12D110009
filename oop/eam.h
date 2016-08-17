#include "datatypes.h"
#include "point3D.h"

#ifndef _EAM_H
  #define _EAM_H

    int* atomicMatrixRead(char* fileName, parameter* p);

    void print_AtomicMatrix(int* mat, int begin, int end);

    Sn_fcc* Sn_fcc_readNeighbours_fromFile(char* fileList);

    void print_Neighbours(Sn_fcc* a);

    Sn_fcc* _defaultFCCNeighbours();

    void eam_data_read(binEAMpot** eam_data, char *fileName, char atom1[2], char atom2[2]);

    rdf* rdf_radius_retrive(binEAMpot* data, double radius);

    eDen_df* eDen_df_charge_density_retrive(binEAMpot* data, double eDen);

    double energyAtIndexFCC(int index, int* a, binEAMpot* data, parameter* p, Sn_fcc* ngbrs);
#endif
