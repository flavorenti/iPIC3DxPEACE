#makefile for PARSEK PROJECT
CPP=	mpicxx
OPTFLAGS=  -O2

SYSTEM = macpro

# UBUNTU 
# include files
INC_HDF5 = -I/usr/local/hdf5/1.8.13/include
# libs
LIB_HDF5 = -L/usr/local/hdf5/1.8.13/lib


#ifeq ($(SYSTEM),supermuc)
# UBUNTU
# include files
#INC_HDF5 = $(HDF5_INC)
# libs
#LIB_HDF5 = ${HDF5_CPP_LIB} ${HDF5_LIB} ${SZIP_LIB}
#endif

#ifeq ($(SYSTEM),ubuntu)
# UBUNTU 
# include files
#INC_HDF5 = -I/home/gianni/computer-stuff/hdf5-1.6.10/hdf5/include
# libs
#LIB_HDF5 = -L/home/gianni/computer-stuff/hdf5-1.6.10/hdf5/lib
#endif

#ifeq ($(SYSTEM),mac)
#MAC
# include files
#INC_HDF5 = -I/Users/gianni/Downloads/hdf5-1.6.10/hdf5/include 
# libs
#LIB_HDF5 = -L/Users/gianni/Downloads/hdf5-1.6.10/hdf5/lib
#endif 

#ifeq ($(SYSTEM),curie)
# CURIE
# include files
#INC_HDF5 = -I/usr/local/hdf5-1.8.8/include
#INC_HDF5 = -I/ccc/cont005/home/ra0747/olshevkv/hdf5-1.6.5/install/include
# libs
#LIB_HDF5 = -L/usr/local/hdf5-1.8.8/lib
#LIB_HDF5 = -L/ccc/cont005/home/ra0747/olshevkv/hdf5-1.6.5/install/lib
#MPELIB = -L/opt/mpi/bullxmpi/1.1.14.3/lib/
#endif

HDF5LIBS =  -lhdf5 -lhdf5_hl  -lz

all: hdf2vtk

hdf2vtk: hdf2vtk_iPic3D_v22mar2016.cpp 
	${CPP} ${OPTFLAGS} -o hdf2vtk ${INC_HDF5} hdf2vtk_iPic3D_v22mar2016.cpp ${LIB_HDF5} ${HDF5LIBS}	
clean:
	rm -rf hdf2vtk



