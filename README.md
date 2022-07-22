@mainpage

# Introduction

### Purpose and goal of the code 
`iPIC3D` was initially conceived by S. Markidis and G. Lapenta in 2010. 
`iPIC3DxPEACE` (iPIC3D for Planets, Exoplanets, Asteroids, Comets and Emissions) is a new code based on `iPIC3D` that targets simulations of the interaction between solar-wind plasma and an obstacle.
Typical case studies are small planetary magnetospheres (like Mercury or Ganymede), exoplanets, magnetized asteroids (like 16 Psyche) and comets (like comet 67P). With this code we can also study emissions (of photons or particles) as a consequence of the plasma-object interaction.

### Structure of the code
`iPIC3DxPEACE` is built upon the state-of-the-art semi-implicit code `iPIC3D` (version KTH 2015) available at this [link](https://github.com/KTH-HPC/iPIC3D).
The semi-implicit numerical scheme is the same as the native `iPIC3D`. Our version uses hybrid parallelization MPI+OpenMP.
The main new features of `iPIC3DxPEACE` are:
1. tailor-made boundary conditions for solar-wind plasma (external) and plasma-planet interaction (internal)
2. time-varying inflow boundary conditions to study the impact of solar-wind fluctuations on magnetospheres (**TBD**)
3. a module to include exospheric neutral atoms and neutral-electron collisions and ionization processes (**ongoing**)
4. a new _in situ_ post-processing module [catalyst](https://www.paraview.org/in-situ/) (**TBD**)
5. extensive documentation generated by [doxygen](https://doxygen.nl) and available at this [link](https://flavorenti.github.io) (**ongoing**)

### Contact information
`iPIC3DxPEACE` is developed and maintained by F. Lavorenti (federico.lavorenti@oca.eu) with contributions from J. Deca and P. Stephenson.
This code was developed at Lagrange laboratory (Observatoire de la Cote d'Azur) and University of Pisa under the PhD project "_Global modelling of Mercury’s outer environment to prepare BepiColombo _" granted by OCA and ESA.
Code available upon request to the authors.

# Installation and User guide

## **Step1** Download the code
`iPIC3DxPEACE` is available upon request to the autors. 
Once you have been granted access to our github [page](https://github.com/flavorenti/iPIC3DxPEACE), you can download the code from a local terminal with the command
```
git clone git@github.com:flavorenti/iPIC3DxPEACE.git
```
This creates a directory `iPIC3DxPEACE` on your computer with a clone copy of the master branch of the code.
Note that since the master branch is secured you can't commit it directly. 
To switch to other branches or create your own new branch you can use the commands `git branch -a` (to show all the branches) and `git branch checkout <name_branch>` (to move to a new or existing branch).

## **Step2** Install dependencies

To compile `iPIC3DxPEACE` you have to have MPI and HDF5 modules on your machine. 
You also need to have a working version of cmake (at least version 2.8), but this can be easily installed on Linux with the command
```
sudo apt-get install cmake
```
For the installation of the dependecies we suggest to follow the User guide developed by CmPA at KU Leuven available at this [link](https://github.com/CmPA/iPic3D/wiki/Quick-User's-Guide). 
To install these modules you don't need sudo or root permissions, everything can be installed locally.
_The H5hut and PETSc modules are not needed for `iPIC3DxPEACE`._

## **Step3** Compile the code

Now that the modules are loaded you can build and run `iPIC3DxPEACE`.
To compile the code run the commands:
```
mkdir build 
cd build 
cmake ..
make
```
If compilation is successful, you will find an executable named **iPIC3D** in build directory. 
This executable is system-dependent but it can be moved around your system to run simulations wherever you want.
To run a test simulation run the commands:
```
cp ../inputfiles/testMagnetosphere3D_small.inp .
mkdir data
mpirun -n 8 ./iPIC3D testMagnetosphere3D_small.inp > RUN.out
```
You will see the text output of the simulation in `RUN.out` and the data in the directory `data`.
If simulation ended successfully: BRAVO! you have run your first `iPIC3DxPEACE` simulation.