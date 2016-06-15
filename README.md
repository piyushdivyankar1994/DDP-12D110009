# DDP-12D110009
Codes for the Dual Degree Project 2016-17 

## 14 June 2016

### Tasks completed till date
1. Neighbour atom coordinates upto 7th nearest neighbours
2. Data Extration code for neighbour atom coordinates
3. EAM data reading function and interpolation function, implemented using binary search.
4. Parameter extraction function
5. Lattice Generator for initial state of montecarlo
6. Lattice reading function
7. Metropolis decision maker

### Incomplete Tasks
1. **Unable to go from crystal coordinate system to lattice index. This makes us unable to compute energy.**
2. EAM energy computation block

## 15 June 2016

### Tasks completed
1. One to one mapping from array to lattice achieved with periodic boundary condition
2. EAM energy calculation done.
3. All monte carlo simulations running without fault
4. Breaking code into functions for better access
5. Simutiion function *eam_monte_carlo_simulation()* is well documented. 
6. New branch created *eam_mc* master reserved for developing other kinds of monteCarlo methods in the future.

### Incomplete Tasks
All assigned tasks completed and commited 

### New Tasks
1. Changing things written as NOTE in the *eam_mc.c* and making the changes
2. Commenting other functions 
3. Write unit tests for all blocks
4. A software to view the fcc atom arrangment data 
5. Device what all anaysis needs to be done and come up with new tasks

### Test Runs
Our code takes in a crystal lattice and EAM data runs montecarlo simulations. We have got this working, so we need to try to generate some results. 

1. Since we are using gsl random number it has a default seed value, and it will always generate the same result. So we need to fix this.
2. We will use parameters Nx=140 Ny=10 Nz=80
3. Lattice parameter 3.56A(Angstrom)
4. We need to write a bash script that will take random seed values generate random input crystals. These input crystal data is subjected to *eam_mc.c* and resulting matrix is stored. 
5. All the input and output data must be saved to seperate folders for future analysis
