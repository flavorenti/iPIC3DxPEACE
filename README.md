/*! @mainpage
 *
 * @section intro_sec Introduction
 *
 * This code was initially conceived by S. Markidis and G. Lapenta in 2010 under the name iPIC3D. 
 * This new version iPIC3DxPEACE (Planets, Exoplanets, Asteroids, Comets and Emissions) targets global simulations of the interaction between the solar wind and an obstacle.
 * iPIC3DxPEACE is developed and maintained by F. Lavorenti with contributions from J. Deca.
 * This code is specifically designed to simulate induced ionospheres and small magnetospheres across the solar system. Specific boundary conditions have been conceived for this purpose.
 * The core of the solver for both particles and fields remains the one of iPIC3D. 
 *
 * @section install_sec Installation
 *
 * @subsection step1 Downloading the code
 *
 * iPIC3DxPEACE is available under-request to the autors. Once you have been granted access to the github page of the code, you can download it by launching the following terminal command
 * @code{.sh}
 * git clone git@github.com:flavorenti/iPIC3D-OCA.git
 * @endcode
 * This command will create a directory iPIC3D-OCA on your computer with a clone copy of the master branch of the code.
 * Note that since the master branch is secured you can't commit it directly. 
 * Instead, you have to create your own branch if you want to modify the code, and then ask for a pull request to merge your own modified branch with the master.
 *
 * @subsection step2 Compiling the code
 *
 * To compile the code you need to install MPI and HDF5 modules. You also need to have a working version of cmake (at least version 2.8), but this can be easily installed with apt-get on Linux systems.
 * Doesn't need sudo or root permissions, everything can be installed locally.
 *
 * 1. Install HDF5 library
 * The current version working properly with iPIC3D is 1.8.9 http://www.hdfgroup.org/ftp/HDF5/releases/
 * List of commands to install HDF5:
 * @code{.sh}
 * cd <top HDF5 source code directory>
 * ./configure --prefix=<location for HDF5 software>
 * make
 * make check
 * make install
 * @endcode
 *
 * 2. Install MPICH2 
 * Install a version that supports -lmpe library. For example, this one: http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/1.3.2p1/mpich2-1.3.2p1.tar.gz
 * During the installation of MPICH follow all the instructions in the README file.
 * Newer versions of OpenMPI should be fine by default in a linux distro.
 *
 * Set environment variables
 * @code{.sh}
 * PATH=mpich2_install_path/bin/:$PATH
 * PATH=hdf5_install_path/lib/:$PATH
 * LD_LIBRARY_PATH=hdf5_install_path/lib/:$LD_LIBRARY_PATH
 * @endcode
 * where mpich2_install_path is the path to the installed MPICH2, and hdf5_install_path is the path to the installed HDF5.
 *
 * @subsection step3 Running the code
 * Create a build directory and move into it.
 * @code{.sh} 
 * cmake ..
 * make
 * @endcode
 * If successful, you will find an executable named iPIC3D in build directory. 
 * This executable is system-dependent but it can be moved around your system to run simulations wherever you want.
 * To run a simulation test case: copy an inputfile named as testXXX.inp from /inputfiles to build directory
 * make sure you create an folder for output as specified in the input file
 * make sure no_of_proc = XLEN x YLEN x ZLEN as specified in the input file
 * @code{.sh}  
 * mpiexec -n no_of_proc ./iPIC3D  inputfilename.inp
 * @endcode
*/

