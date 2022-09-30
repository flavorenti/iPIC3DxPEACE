#!/bin/bash

mkdir build

cd build

cmake ..

make

cp ../inputfiles/*_small.inp . 

mkdir data

mpirun -n 4 ./iPIC3D testEmptyBox3D_small.inp

echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30

mpirun -n 8 ./iPIC3D testMagnetosphere3D_small.inp
