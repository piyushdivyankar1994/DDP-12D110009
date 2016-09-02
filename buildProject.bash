 #!

if gcc -std=c99 -Wall -I./include main.c src/simulation.c src/chemicalPotential.c src/analysis.c src/eam.c src/parameters.c src/math_functions.c src/point3D.c src/test.c -o test.out -lm -lgsl -lgslcblas -g; then

    echo "######################################################################"
    echo "Program compliled now running valgrind"
    echo "######################################################################"
    valgrind ./test.out --leak-check=full --show-reachable=no --track-origins=yes  --trace-children=yes -v --num-callers=100
    #echo "do you want to run gdb (y/n)"
fi
