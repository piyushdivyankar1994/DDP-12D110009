#include<stdio.h>
#include<stdlib.h>
#include "math_functions.h"
#include "string.h"

/** TO implement random numbers*/
#include "math.h"
#include "gsl/gsl_rng.h"

gsl_rng * r;


#ifndef _DATATYPES_H
#define _DATATYPES_H 1
    struct parameter {
      int no_of_atoms;
      int Nx;
      int Ny;
      int Nz;
      double lattice_parameter;
      int N_MCS;
      double temperature;
      int atoms_per_site;
      int nearestNeighbours[7];
    }typedef parameter;

    struct binaryEAMpotential {
      char atom1[2];
      int no_of_files;
      char atom2[2];
      double radius[5000];
      double atom1_charge_Density[5000];
      double atom2_charge_Density[5000];
      double pair_atom1_atom1[5000];
      double pair_atom1_atom2[5000];
      double pair_atom2_atom2[5000];
      double atom1_embedding_energy[5000];
      double atom2_embedding_energy[5000];
      double chargeDensity[5000];
      double minRadius;
      double maxRadius;
      double min_eDen;
      double max_eDen;
    }typedef binEAMpot;

    struct radius_dependent_fields {
      int index;
      double radius;
      double p11;
      double p12;
      double p22;
      double eDen1;
      double eDen2;
    }typedef rdf;

    struct chargeDensity_dependent_fields {
      int index;
      double eDen;
      double embed1;
      double embed2;
    }typedef eDen_df;

    struct point3D {
      float x;
      float y;
      float z;
    }typedef point3D;


    struct neighbours_fcc{
    	point3D s1n[12];
    	point3D s2n[6];
    	point3D s3n[24];
    	point3D s4n[12];
    	point3D s5n[24];
    	point3D s6n[8];
    	point3D s7n[48];
    	int indices[7];
    }typedef Sn_fcc;

    enum Boolean {
      true,
      false
    }typedef Boolean;

#endif
