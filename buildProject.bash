 #!

gcc -std=c99 -Wall -I./include main.c src/simulation.c src/chemicalPotential.c src/analysis.c src/eam.c src/parameters.c src/math_functions.c src/point3D.c src/test.c -o test.out -lm -lgsl -lgslcblas -g

./test.out
