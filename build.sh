#!/bin/bash


echo "STARTED BUILDING C++ COMPONENTS"
echo "USING COMPILER"
echo `which g++`
echo `g++-9 -v`

[ -d bin ] && rm -r bin
mkdir bin

echo "BUILDING THE PlasLR 15 MER COMPUTATIONS"
g++-9 src/search-15mers.cpp -fopenmp -Wall -o bin/search15mers
echo "BUILDING THE PlasLR 3 MER COMPUTATIONS"
g++-9 src/count-tri.cpp -Wall -fopenmp -o bin/countTrimers
echo "BUILDING THE PlasLR READ FILTER"
g++-9 src/filter_reads.cpp -Wall -fopenmp -o bin/filter

echo "BUILD FINISHED"
