cmake_minimum_required(VERSION 2.8.8) 

#
# Project declaration
#
project(iPic3D)

#
# Set compiler flags per system
#
# Possible flags are the following:
# For Xeon Phi, totalview: -openmp -g -mmic
# Reporting: -g -vec-report
# Optimization: -O3 -xHost -fno-exceptions
# OpenMP: -openmp
# Print pcls subcycles: -DPRINTPCL 
# Disable HDF5 output: -DNO_HDF5
set(CMAKE_CXX_FLAGS "-O3")

#
# Set on which computer we are running
#
# possible values are IRENE, LICALLO, LOCAL
set(USE_MACHINE "LOCAL")

#
# Catalyst
#
# USE_CATALYST enables catalyst. Possible values are
# "LEGACY" or "V2". Every other value will turn it off.
set(USE_CATALYST "NO")

#
# Catalyst Legacy - load modules and libraries
#
if(USE_CATALYST STREQUAL "LEGACY")

  add_definitions(-DUSE_CATALYST_LEGACY)    

  if(USE_MACHINE STREQUAL "LOCAL")
    set(ParaView_DIR /home/flavorenti/Desktop/Sources/paraview_lapo/build) #path to paraview build
  endif()
  if(USE_MACHINE STREQUAL "LICALLO")
    set(ParaView_DIR /data/home/lquerci/Sources/paraview/build)     
  endif()
  if(USE_MACHINE STREQUAL "IRENE")
    list(APPEND CMAKE_MODULE_PATH "/ccc/products/paraview-5.8.0/intel--19.0.5.281__openmpi--4.0.1/default/lib64/cmake/paraview-5.8")
    set(ParaView_DIR /ccc/products/paraview-5.8.0/intel--19.0.5.281__openmpi--4.0.1/default/lib64/cmake/paraview-5.8)
  endif()

  find_package(ParaView 5.7 REQUIRED)

  if (NOT TARGET ParaView::PythonCatalyst)
    message(STATUS
      "${CMAKE_PROJECT_NAME} requires ParaView to be built with Catalyst and \n"
      "Python support enabled. Please rebuild ParaView (or point to a \n"
      "different build of ParaView) with PARAVIEW_ENABLE_CATALYST and \n"
      "PARAVIEW_ENABLE_PYTHON set to TRUE")
  else()
    add_library(iPICAdaptor ${CMAKE_SOURCE_DIR}/catalyst_legacy/Adaptor_legacy.cpp)
    if(USE_MACHINE STREQUAL "IRENE")
      target_link_libraries(iPICAdaptor PRIVATE ParaView::PythonCatalyst VTK::CommonDataModel Python2::Python)
    else()
      target_link_libraries(iPICAdaptor PRIVATE ParaView::PythonCatalyst VTK::CommonDataModel Python3::Python)
    endif()
    target_include_directories(iPICAdaptor INTERFACE ${CMAKE_SOURCE_DIR}/catalyst_legacy)
    target_compile_definitions(iPICAdaptor INTERFACE USE_CATALYST)
  endif()
  
endif()

#
# Find MPI library
#
if(USE_MACHINE STREQUAL "IRENE")
   cmake_policy(SET CMP0074 NEW)
endif()
find_package(MPI)
set(EXTRA_LIBS "")   

#
# Find HDF5 library
#
if(${CMAKE_CXX_FLAGS} MATCHES "(.*)-DNO_HDF5(.*)")
# do nothing
else()
  if(USE_MACHINE STREQUAL "LICALLO")
    cmake_policy(SET CMP0074 NEW)
  endif()
  add_definitions(-DPHDF5)
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 COMPONENTS HL C )   
endif()

#
# include and lib directories
#
include_directories(
        include
        ${HDF5_INCLUDE_DIRS}
        ${HDF5_HL_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

#
# Link libraries - if machine=Irene need to explicitly link intel63 compiler
# 
if(USE_MACHINE STREQUAL "IRENE")
  link_directories(
      ${HDF5_LIBRARY_DIRS}
      ${HDF5_hdf5_hl_LIBRARY_DIRS}
      ${MPI_LIBRARY_DIRS}
      ${EXTRA_LIBS}
      /ccc/products/icc-20.0.4/system/default/20.0.4/lib/intel64
  )
else()
  link_directories(
      ${HDF5_LIBRARY_DIRS}
      ${HDF5_hdf5_hl_LIBRARY_DIRS}
      ${MPI_LIBRARY_DIRS}
      ${EXTRA_LIBS}  
  )
endif()

#
# Header file list
#
file(
        GLOB
        inc_files
        include/*.h
)

#
# Source file list
#
file(
        GLOB
        src_files
        ConfigFile/*.cpp
        PSKOutput3D/*.cpp
        bc/*.cpp
        communication/*.cpp
        fields/*.cpp
        grids/*.cpp
        inputoutput/*.cpp
        particles/*.cpp
        performances/*.cpp
        solvers/*.cpp
        utility/*.cpp
        main/*.cpp
)

#
# Executable declaration
#
add_executable(
        iPIC3D
        iPIC3D.cpp
)
option(STATIC_LINK "Choose static link"  on)
if(STATIC_LINK) 
	add_library(
	    iPIC3Dlib       # name of the library
		  STATIC    # type of the library
	    ${inc_files}    # stuff to build the library
	    ${src_files}
	)
else()
	add_library(
	    iPIC3Dlib       # name of the library
		  SHARED    # type of the library
	    ${inc_files}    # stuff to build the library
	    ${src_files}
	)
endif()


#
# Link external libraries
#
if (USE_MACHINE STREQUAL "IRENE")
  if (USE_CATALYST STREQUAL "LEGACY")
    target_link_libraries(
	 iPIC3Dlib
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
	 ${EXTRA_LIBS}
	 /ccc/products/icc-20.0.4/system/default/20.0.4/lib/intel64/libiomp5.so
    	 /ccc/products2/icc-19.0.5.281/Atos_7__x86_64/system/default/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin/libirc.so
  )
  else()
    target_link_libraries(
	 iPIC3Dlib
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
	 ${EXTRA_LIBS}
	 /ccc/products/icc-20.0.4/system/default/20.0.4/lib/intel64/libiomp5.so
    )
  endif()
else()
  target_link_libraries(
	iPIC3Dlib
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
	 ${EXTRA_LIBS}
  )
endif()


target_link_libraries(
         iPIC3D
         iPIC3Dlib
)

if (USE_CATALYST STREQUAL "LEGACY")
  target_link_libraries(
      iPIC3Dlib
      iPICAdaptor
  )
endif()

#
# Properties of executable
#
# to save the executable in the folder 
# where the CMakeLists.txt file is use CMAKE_CURRENT_SOURCE_DIR
# debug releases have a _d appended to the executable
set_target_properties(iPIC3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(iPIC3D PROPERTIES DEBUG_POSTFIX "_d")

#
# Final message at the end of compilation
#
message("Which system am I compiling for:")
message("USER: $ENV{USER}")
message("HOSTNAME: $ENV{HOSTNAME}")
message("CMAKE_SYSTEM_PROCESSOR is ${CMAKE_SYSTEM_PROCESSOR}")
message("Compiler & compiler flags:")
message("CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("HDF5_INCLUDE_DIRS is ${HDF5_INCLUDE_DIRS}")
message("HDF5_LIBRARIES is ${HDF5_LIBRARIES}")
message("HDF5_HL_LIBRARIES is ${HDF5_HL_LIBRARIES}")
message("MPI_LIBRARIES is ${MPI_LIBRARIES}")
message("USE_CATALYST is set to ${USE_CATALYST}")
message("USE_MACHINE is set to ${USE_MACHINE}")
