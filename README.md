# Monte-Carlo simulation of interfacial energy calculations in soild-solid systems
(Dual Degree Project - IIT Bombay)

### Author: Piyush Divyankar
### Guide: Prof. M.P. Gururajan

## Working functions
- Embedded atom method can be used to calculate energies in FCC type fixed lattices. The examples of this can be seen in **src/simulation.c**.
- Fixed lattice simulations possible for FCC type crystals.
- Neighbour index table, i.e a table that maps each index to indexes of neighbouring points.
- EAM energy table that stores energies of all possible configurations around a point.
- A small point3D library for perfoming functions related to miller indices.

## To Do

- Make more functions in analysis.c as need arises.
- Generate relevant physical results.

## Objectives completed.
- Code benchmarked for Ni-Al system in lattice parameter.
- Ni-Al system cannoical ensemble benchmarked using LJP type potential.
- Required functionality for dual degree project is achieved.

## Future Goals
- New brach created general systems.
- The code should be modified to fit any crystal that can be broken to a equivalent rectangular coordiante system(except hexagonal) with any motiff.
- All functionality available for FCC must be applicable to general systems. This inculdes but not limited to
  - Neighbour Index table.
  - Energy table for EAM
  - Periodic boundary transform
