# Details about some stuff

## Neighbour sites files

These files are named S1n.mat, S2n.mat etc. The number in the middle of filename indicates n'th nearest neighbour. Deleting these will cause code to break down.

## EAM data files 

Names of all EAM data files is maintained in file *input_data.txt*. Names of these files are like fAl.. F_Ni.. pAl.. . Absolutely essential information

## Parameters file

Parameters file must have first entry an integer. This will maintain total no of useful values in the file. The next n lines must contain information. Except for the first rest all can be real numbers. Any additional data will be lost. If there as less values then the array created in the program will contain garbage values.
If above convention isn't followed to the letter can lead to segmentaion fault. 

### Current *parameters.txt*
It contains following information
1. No of entries
2. Nx
3. Ny
4. Nz
5. Disordered phase percentage
6. Ordered phase percentage
7. Percentage of nickel in disordered phase. 

## FCC crystal input data
A little about how we store the fcc crystal 

Every FCC lattice can be thought of a simple cubic lattice with a four atom motiff.
Motiff atoms have following positions .

(0.0 0.0 0.0),(0.0 0.5 0.5),(0.5 0.0 0.5),(0.5 0.5 0.0)

We store the species information in a linear array such that atoms associated with each lattice points are stored sucessively 
the motiff atoms are needs to acessed uniquely so a coressponds to which atom we are taking about
a=0 , 1 , 2 , 3 represent the 4 types of location that an atom can have

.this data is stored in *input_data.txt*
 
output data format is same in the file

## Random Initial Step 

By this we mean, *matrix_generator.c* can generate a random lattice that will serve as an random initial step for the Montecarlo process. This function can run using the bash script *matrix_compile_and_run_once*
