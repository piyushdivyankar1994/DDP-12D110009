/*
   This file is part of DDP-12D110009.

    DDP-12D110009 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DDP-12D110009 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

/*!
   \file datatypes.h
   \brief This file contains all composite datatypes used in this Project
   \author Piyush Divyankar
   \date 30/08/2016
 */
#include <stdio.h>
#include <stdlib.h>
#include "math_functions.h"
#include "string.h"

/** TO implement random numbers*/
#include "math.h"
#include "gsl/gsl_rng.h"

gsl_rng * r;
/**
 * Default name for parameter file.
 */
#define PARAM_FILE_DEF_NAME "defaultFCC.param"

#ifndef _DATATYPES_H
#define _DATATYPES_H        1

/**
 * This structure contains list of parameters that define the system.
 * \struct parameter
 * \typedef parameter
 * \author Piyush Divyankar
 */
struct parameter
{
    /** Total number of atoms in the simulation matrix */
    int no_of_atoms;
    /** Size of lattice in x direction */
    int Nx;
    /** Size of lattice in y direction */
    int Ny;
    /** Size of lattice in z direction */
    int Nz;
    /** Size of a unit cell in the lattice */
    double lattice_parameter;
    /** Total number of Monte-Carlo simulations per atom */
    int N_MCS;
    /** temperature of material at time of simulation */
    double temperature;
    /** Number of atoms per site in motiff */
    int atoms_per_site;
    /** Number of nth nearest neighbours stores upto 7 */
    // / FUTURE_CHANGES:20 must make this general so that upto nth nearest neighbours can be
    // / shown
    int nearestNeighbours[7];
    /** Stores the FILE name that the object is parsed from unset from default*/
    /** TODO:0 make so that it is '\0' by default */
    /** DONE:20 change functions of parsing this type to reflect this addition */
    char fileName[50];
} typedef parameter;
/**
 * This structure is meant to store EAM data for energy computation.
 * \struct binaryEAMpotential
 * \typedef binEAMpot
 * \author Piyush Divyankar
 * \note This structure currently conforms to data that is placed in /EAM_Ni_Al/
 *       In future changes must be made so that this size isn't fixed.
 */
struct binaryEAMpotential
{
    /** Chemical symbol for atom 1 involved in simulation */
    char atom1[2];
    /** Chemical symbol for atom 2 in simulation */
    char atom2[2];
    /** Total no of EAM data files */
    int no_of_files;
    /**
     * Key value radius. This is used to index other values that are a function
     * of radius. This remains same across various files, that contain radius
     * dependent data. This was the case for EAM data that I am using for proof
     * of concept. It might be different later.
     */
    double radius[5000];
    /** Charge density due to atom 1 at given distance from atom */
    double atom1_charge_Density[5000];
    /** Charge density due to atom 2 at given distance from atom */
    double atom2_charge_Density[5000];
    /** Electrostatic pair wise potential for atom1 at given distance */
    double pair_atom1_atom1[5000];
    /** Electrostatic pair wise potential for atom1 and atom2 at given distnace */
    double pair_atom1_atom2[5000];
    /** Electrostatic pair wise potential for atom2 at given distance */
    double pair_atom2_atom2[5000];
    /**
     * Key value charge density. Embedding energies are functions of charge
     * density at a given location, and not radius. This array holds various
     * values at which embedding energies are calculated.
     */
    double chargeDensity[5000];
    /** Embedding energy for atom 1 at given charge density */
    double atom1_embedding_energy[5000];
    /** Embedding energy for atom 2 at given charge density */
    double atom2_embedding_energy[5000];
    /** Minimum value of radius usually the 0 index */
    double minRadius;
    /** Maximum value of radius usually the 5000 index */
    double maxRadius;
    /** Minimum value of charge density usually the 0 index */
    double min_eDen;
    /** Maximum value of charge density usually the 5000 index */
    double max_eDen;
} typedef binEAMpot;

/**
 * This structure is meant to store EAM data of a single row of values that
 * depend on radius for energy computation.
 * \struct radius_dependent_fields
 * \typedef rdf
 * \author Piyush Divyankar
 */
struct radius_dependent_fields
{
    /** Index of the record that was picked from the table */
    int index;
    /** Key value of the radius */
    double radius;
    /** Linearly interpolated pairwise potential atom1 */
    double p11;
    /** Linearly interpolated pairwise potential atom1-atom2 */
    double p12;
    /** Linearly interpolated pairwise potential atom2 */
    double p22;
    /** Linearly interpolated charge density atom 1 */
    double eDen1;
    /** Linearly interpolated charge density atom 2 */
    double eDen2;
} typedef rdf;

/**
 * This structure is meant to store EAM data of a single row of values that
 * depend on charge density for energy computation.
 * \struct radius_dependent_fields
 * \typedef eDen_df
 * \author Piyush Divyankar
 */
struct chargeDensity_dependent_fields
{
    /** Index of the values just smaller than key value */
    int index;
    /** Key value */
    double eDen;
    /** Linearly interpolated embedding energy of atom 1 at key value */
    double embed1;
    /** Linearly interpolated embedding energy of atom 2 at key value */
    double embed2;
} typedef eDen_df;

/**
 * This structure is used to represent 3D millier indices of points in real space
 * of the crystal lattice.
 * \struct point3D
 * \typedef point3D
 * \author Piyush Divyankar
 */
struct point3D
{
    /** x-coordinate */
    float x;
    /** y-coordinate */
    float y;
    /** z-coordinate */
    float z;
} typedef point3D;

/**
 * This structure holds an array of points that represent the relative miller
 * indices of the neighbours atoms from a particular site.
 * \struct neighbours_fcc
 * \typedef Sn_fcc
 * \author Piyush Divyankar
 */
struct neighbours_fcc
{
    /** List of points for first nearest neighbours */
    point3D s1n[12];
    /** List of points for second nearest neighbours */
    point3D s2n[6];
    /** List of points for third nearest neighbours */
    point3D s3n[24];
    /** List of points for fourth nearest neighbours */
    point3D s4n[12];
    /** List of points for fifth nearest neighbours */
    point3D s5n[24];
    /** List of points for sixth nearest neighbours */
    point3D s6n[8];
    /** List of points for seventh nearest neighbours */
    point3D s7n[48];
    /** Contains the numbers of first neareast, second neareast etc. */
    int indices[7];
} typedef Sn_fcc;

/** \note Not really sure why I did this. Could be useful in some ways */
enum Boolean
{
    true,
    false
} typedef Boolean;
/**
 * \typedef ATOM
 * \brief Atoms are stored as integers in memory. So, just to avoid confusion.
 */
typedef int ATOM;

#endif /* ifndef _DATATYPES_H */
