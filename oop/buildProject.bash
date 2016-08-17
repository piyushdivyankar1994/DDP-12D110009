#!

gcc -std=c99 -Wall main.c chemicalPotential.c eam.c parameters.c math_functions.c point3D.c test.c -o test.out -lm -g

./test.out
