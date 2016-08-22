#!

gcc -std=c99 -Wall main.c simulation.c chemicalPotential.c analysis.c eam.c parameters.c math_functions.c point3D.c test.c -o test.out -lm -lgsl -lgslcblas -g

./test.out
