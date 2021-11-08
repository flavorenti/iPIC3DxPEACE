#!/bin/bash

mkdir build

cd build

cmake ..

make

cp ../inputfiles/testMagnetosphere3Dsmall.inp ../inputfiles/testMagnetosphere2Dsmall.inp .

mkdir data

mpirun -n 8 ./iPIC3D testMagnetosphere3Dsmall.inp
mpirun -n 6 ./iPIC3D testMagnetosphere2Dsmall.inp
