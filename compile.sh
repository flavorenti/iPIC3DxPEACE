#!/bin/bash

mkdir build

cd build

cmake ..

make

cp ../inputfiles/testMagnetosphere3D_small.inp ../inputfiles/testMagnetosphere2D_small.inp .

mkdir data

mpirun -n 8 ./iPIC3D testMagnetosphere3D_small.inp
mpirun -n 4 ./iPIC3D testMagnetosphere2D_small.inp
