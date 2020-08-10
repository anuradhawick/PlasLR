#!/bin/bash


echo "STARTED BUILDING C++ COMPONENTS"
echo "USING COMPILER"
echo `which g++`
echo `g++ -v`

[ -d bin ] && rm -r bin
mkdir bin

echo "BUILDING THE PlasLR 15 MER COMPUTATIONS"
g++ src/search-15mers.cpp -fopenmp -Wall -o bin/search15mers
echo "BUILDING THE PlasLR 3 MER COMPUTATIONS"
g++ src/count-tri.cpp -Wall -fopenmp -o bin/countTrimers

echo "BUILD FINISHED"
