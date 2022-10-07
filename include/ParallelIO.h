/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef __PARALLELIO_H__
#define __PARALLELIO_H__

#ifdef USEH5HUT
#  include "H5hut-io.h"
#endif

#ifdef PHDF5
#  include "phdf5.h"
#endif

#include "ipicfwd.h"
#include "arraysfwd.h"
#include <string>
using std::string;

void WriteFieldsH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle);
void WritePartclH5hut(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);

void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid);
void ReadFieldsH5hut(int nspec, EMfields3D *EMf,       Collective *col, VCtopology3D *vct, Grid3DCU *grid);

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);



/**************MPI_IO*********************/
int WriteFieldsVTKNonblk(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct,int cycle,
			float**** fieldwritebuffer,MPI_Request requestArr[4],MPI_File fhArr[4]);

int WriteMomentsVTKNonblk(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct,int cycle,
			float*** momentswritebuffer,MPI_Request requestArr[14],MPI_File fhArr[14]);

void WriteFieldsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle);
void WriteFieldsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,float**** fieldwritebuffer);
void WriteMomentsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,float***  momentswritebuffer);

/*! @brief Write the energy spectra into vtk files.
 *
 * The function takes the volume of the cells and the desired energy bins from the input configuration file
 * It cycles through all the particles of the simulation, finds the cell where each particle is and it sums its 
 * charge to the corresponding energy bin, defining the energy spectra on each cell.
 * It defines the energy spectra with velocities parallel to the magnetic field, perpendicular and total.
 * It creates vtk files with the the absolute total charge of the species in each energy bin of each cell.
 *
 * @param[in] grid built-in type defining the uniform cartesian 3D local grid for each process
 * @param[in] part built-in type for particles of the same species, in a 3D space and 3 component velocity
 * @param[in] EMf built-in type with the electromagnetic fields and sources defined for each local grid
 * @param[in] col built-in type with all the collective properties of the simulation
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] tag tags for the desired files from the input configuration file (Spar, Sperp, Stot)
 * @param[in] cycle number of the current iteration
 * @param[in] spectrawritebuffere 4D array for electrons: [x index of cell][y index of cell][z index of cell][total absolute charge in different bins on the cell]
 * @param[in] spectrawritebufferi 4D array for ions: [x index of cell][y index of cell][z index of cell][total charge in different bins on the cell]
 */
void WriteSpectraVTK(Grid3DCU *grid, Particles3D *part, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,float****  spectrawritebuffere,float****  spectrawritebufferi);

/*! @brief Write the temperature tensor into vtk files.
 *
 * The function calculates the temperature tensor in cartesian coordinates and coordinates
 * relative to the magnetic field where the diagonals of the pressure tensor are in gyrotropic form.
 * It creates vtk files with the 6 components of the symmetric tensor in all the cells.
 *
 * @param[in] grid built-in type defining the uniform cartesian 3D local grid for each process
 * @param[in] EMf built-in type with the electromagnetic fields and sources defined for each local grid
 * @param[in] col built-in type with all the collective properties of the simulation
 * @param[in] vct built-in type defining the 3D MPI cartesian topology
 * @param[in] tag tags for the desired files from the input configuration file (Tcart, Tperpar)
 * @param[in] cycle number of the current iteration
 * @param[in] temperaturewritebuffer 4D array: [x index of cell][y index of cell][z index of cell][components of the tensor]
 */
void WriteTemperatureVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,float**** temperaturewritebuffer);
void WriteTestPclsVTK(int nspec, Grid3DCU *grid, Particles3D *part, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle, MPI_Request *testpartMPIReq, MPI_File *fh);
void ByteSwap(unsigned char * b, int n);
#endif
