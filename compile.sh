#!/bin/bash

mkdir build

cd build

cmake ..

make

cp ../inputfiles/*_small*.inp . 

mkdir data

mpirun -n 8 ./iPIC3D testEmptyBox3D_small_Pete.inp

echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30

mv data data1
mkdir data

mpirun -n 8 ./iPIC3D testMagnetosphere3D_small.inp

echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30

mv data data2
mkdir data

mpirun -n 8 ./iPIC3D testMagnetosphere3D_small_Julen.inp

mv data data3
