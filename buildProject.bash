 #!

if gcc -std=c99 -Wall -I./include main.c src/simulation.c src/chemicalPotential.c src/analysis.c src/eam.c src/parameters.c src/math_functions.c src/point3D.c src/test.c -o test1.out -lm -lgsl -lgslcblas -g -pg; then

    echo "######################################################################"
    echo "Program compliled"
    echo "######################################################################"
    #valgrind ./test.out --leak-check=full --show-reachable=no --track-origins=yes  --trace-children=yes -v --num-callers=100
    #echo "do you want to run gdb (y/n)"

fi
echo 'Program Running'
echo "######################################################################"

./test1.out > log.txt

echo 'Plotting'
echo "######################################################################"

echo 'load "plotfile"' | gnuplot

xdg-open lattice_parameter_convergence.png
cat log.txt
