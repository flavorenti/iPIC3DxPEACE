#!/bin/bash

# Compile the code into build
mkdir build
cd build
cmake ..
make

# Copy test inputfiles to build
cp ../inputfiles/*_small*.inp . 

# 1st test with empty box and no collision. Typical SW conditions found at Mercury.
mkdir data
mpirun -n 8 ./iPIC3D testEmptyBox3D_noColl_small.inp > RUN_data1.out
mv data/ConservedQuantities.txt data/ConservedQuantities0.txt
mpirun -n 8 ./iPIC3D restart testEmptyBox3D_noColl_small.inp > RUN_data1.res
mv data/ConservedQuantities.txt data/ConservedQuantities1.txt
echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30
mv data data1

# 2nd test with empty box and yes collision. Uniform neutral density in the box.
mkdir data
mpirun -n 8 ./iPIC3D testEmptyBox3D_yesColl_small.inp > RUN_data2.out
mv data/ConservedQuantities.txt data/ConservedQuantities0.txt
mpirun -n 8 ./iPIC3D restart testEmptyBox3D_yesColl_small.inp > RUN_data2.res
mv data/ConservedQuantities.txt data/ConservedQuantities1.txt
echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30
mv data data2

# 3rd test with collisions on a localized cloud of neutrals.
mkdir data
mpirun -n 8 ./iPIC3D testExosphere3D_yesColl_small.inp > RUN_data3.out
mv data/ConservedQuantities.txt data/ConservedQuantities0.txt
mpirun -n 8 ./iPIC3D restart testExosphere3D_yesColl_small.inp > RUN_data3.res
mv data/ConservedQuantities.txt data/ConservedQuantities1.txt
echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30
mv data data3

# 4th test with ionization of the cloud of neutrals and no collisions.
mkdir data
mpirun -n 8 ./iPIC3D testExosphere3D_yesIoniz_small.inp > RUN_data4.out
mv data/ConservedQuantities.txt data/ConservedQuantities0.txt
mpirun -n 8 ./iPIC3D restart testExosphere3D_yesIoniz_small.inp > RUN_data4.res
mv data/ConservedQuantities.txt data/ConservedQuantities1.txt
echo "LETS WAIT 30 SEC BEFORE RUNNING THE NEW TEST..."
sleep 30
mv data data4

# 5th test with both ionization and collisions on the cloud of neutrals.
mkdir data
mpirun -n 8 ./iPIC3D testExosphere3D_yesAll_small.inp > RUN_data5.out
mv data/ConservedQuantities.txt data/ConservedQuantities0.txt
mpirun -n 8 ./iPIC3D restart testExosphere3D_yesAll_small.inp > RUN_data5.res
mv data/ConservedQuantities.txt data/ConservedQuantities1.txt
mv data data5

# remove the inputfiles from the build directory
rm ./*_small*.inp
