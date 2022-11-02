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

#include <mpi.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "ParallelIO.h"
#include "MPIdata.h"
#include "TimeTasks.h"
#include "Collective.h"
#include "Grid3DCU.h"
#include "VCtopology3D.h"
#include "Particles3D.h"
#include "EMfields3D.h"
#include "math.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>

/*! Function used to write the EM fields using the parallel HDF5 library */
void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef PHDF5
  timeTasks_set_task(TimeTasks::WRITE_FIELDS);

  stringstream filenmbr;
  string       filename;

  bool         bp;

  /* ------------------- */
  /* Setup the file name */
  /* ------------------- */

  filenmbr << setfill('0') << setw(5) << cycle;
  filename = col->getSaveDirName() + "/" + col->getSimName() + "_" + filenmbr.str() + ".h5";

  /* ---------------------------------------------------------------------------- */
  /* Define the number of cells in the globa and local mesh and set the mesh size */
  /* ---------------------------------------------------------------------------- */

  int nxc = grid->getNXC();
  int nyc = grid->getNYC();
  int nzc = grid->getNZC();

  int    dglob[3] = { col ->getNxc()  , col ->getNyc()  , col ->getNzc()   };
  int    dlocl[3] = { nxc-2,            nyc-2,            nzc-2 };
  double L    [3] = { col ->getLx ()  , col ->getLy ()  , col ->getLz ()   };

  /* --------------------------------------- */
  /* Declare and open the parallel HDF5 file */
  /* --------------------------------------- */

  PHDF5fileClass outputfile(filename, 3, vct->getCoordinates(), vct->getFieldComm());

  bp = false;
  if (col->getParticlesOutputCycle() > 0) bp = true;

  outputfile.CreatePHDF5file(L, dglob, dlocl, bp);

  // write electromagnetic field
  //
  outputfile.WritePHDF5dataset("Fields", "Ex", EMf->getEx(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Ey", EMf->getEy(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Ez", EMf->getEz(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Bx", EMf->getBx(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "By", EMf->getBy(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Bz", EMf->getBz(), nxc-2, nyc-2, nzc-2);

  /* ---------------------------------------- */
  /* Write the charge moments for each species */
  /* ---------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream ss;
    ss << is;
    string s_is = ss.str();

    // charge density
    outputfile.WritePHDF5dataset_mom("Fields", "rho_"+s_is, EMf->getRHOcs(), is, nxc-2, nyc-2, nzc-2);
    // current
    outputfile.WritePHDF5dataset_mom("Fields", "Jx_"+s_is, EMf->getJxs(), is, nxc-2, nyc-2, nzc-2);
    outputfile.WritePHDF5dataset_mom("Fields", "Jy_"+s_is, EMf->getJys(), is, nxc-2, nyc-2, nzc-2);
    outputfile.WritePHDF5dataset_mom("Fields", "Jz_"+s_is, EMf->getJzs(), is, nxc-2, nyc-2, nzc-2);
  }

  outputfile.ClosePHDF5file();

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}

/*! Function to write the EM fields using the H5hut library. */
void WriteFieldsH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){
  if(col->field_output_is_off())
    return;
#ifdef USEH5HUT
  timeTasks_set_task(TimeTasks::WRITE_FIELDS);

  H5output file;


  /* ---------------- */
  /* Write the fields */
  /* ---------------- */

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetNameCycle(filename, cycle);

  file.OpenFieldsFile("Node", nspec, col->getNxc()+1, col->getNyc()+1, col->getNzc()+1, vct->getCoordinates(), vct->getDims(), vct->getFieldComm());

  file.WriteFields(EMf->getEx(), "Ex", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getEy(), "Ey", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getEz(), "Ez", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBx(), "Bx", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBy(), "By", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBz(), "Bz", grid->getNXN(), grid->getNYN(), grid->getNZN());

  for (int is=0; is<nspec; is++) {
    stringstream  ss;
    ss << is;
    string s_is = ss.str();
    file.WriteFields(EMf->getRHOns(is), "rho_"+ s_is, grid->getNXN(), grid->getNYN(), grid->getNZN());
  }

  file.CloseFieldsFile();

//--- SAVE FIELDS IN THE CELLS:
//
//  file.OpenFieldsFile("Cell", nspec, col->getNxc(), col->getNyc(), col->getNzc(), vct->getCoordinates(), vct->getDims(), vct->getComm());
//
//  file.WriteFields(EMf->getExc(), "Exc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getEyc(), "Eyc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getEzc(), "Ezc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getBxc(), "Bxc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getByc(), "Byc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getBzc(), "Bzc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//
//  for (int is=0; is<nspec; is++) {
//    stringstream  ss;
//    ss << is;
//    string s_is = ss.str();
//    file.WriteFields(EMf->getRHOcs(is), "rhoc_"+ s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
//  }
//
//  file.CloseFieldsFile();
//
//--- END SAVE FIELDS IN THE CELLS.

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}

/*! Function to write the particles using the H5hut library. */
void WritePartclH5hut(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){
#ifdef USEH5HUT
  timeTasks_set_task(TimeTasks::WRITE_PARTICLES);

  H5output file;

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetNameCycle(filename, cycle);

  /* ------------------- */
  /* Write the particles */
  /* ------------------- */

  file.OpenPartclFile(nspec, vct->getFieldComm());
  for (int i=0; i<nspec; i++){
    // this is a hack
    part[i].convertParticlesToSynched();
    file.WriteParticles(i, part[i].getNOP(),
                           part[i].getQall(),
                           part[i].getXall(),
                           part[i].getYall(),
                           part[i].getZall(),
                           part[i].getUall(),
                           part[i].getVall(),
                           part[i].getWall(),
                           vct->getFieldComm());
  }
  file.ClosePartclFile();

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}

#if 0
void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid){
#ifdef USEH5HUT

  H5input infile;
  double L[3] = {col->getLx(), col->getLy(), col->getLz()};

  infile.SetNameCycle(col->getinitfile(), col->getLast_cycle());
  infile.OpenPartclFile(nspec);

  infile.ReadParticles(vct->getCartesian_rank(), vct->getNproc(), vct->getDims(), L, vct->getFieldComm());

  for (int s = 0; s < nspec; s++){
    part[s].allocate(s, infile.GetNp(s), col, vct, grid);

    infile.DumpPartclX(part[s].getXref(), s);
    infile.DumpPartclY(part[s].getYref(), s);
    infile.DumpPartclZ(part[s].getZref(), s);
    infile.DumpPartclU(part[s].getUref(), s);
    infile.DumpPartclV(part[s].getVref(), s);
    infile.DumpPartclW(part[s].getWref(), s);
    infile.DumpPartclQ(part[s].getQref(), s);
  }
  infile.ClosePartclFile();

//--- TEST PARTICLE LECTURE:
//  for (int s = 0; s < nspec; s++){
//    for (int n = 0; n < part[s].getNOP(); n++){
//      double ix = part[s].getX(n);
//      double iy = part[s].getY(n);
//      double iz = part[s].getZ(n);
//      if (ix<=0 || iy<=0 || iz <=0) {
//        cout << " ERROR: This particle has negative position. " << endl;
//        cout << "        n = " << n << "/" << part[s].getNOP();
//        cout << "       ix = " << ix;
//        cout << "       iy = " << iy;
//        cout << "       iz = " << iz;
//      }
//    }
//  }
//--- END TEST

#endif
}
#endif

#if 0
void ReadFieldsH5hut(int nspec, EMfields3D *EMf, Collective *col, VCtopology3D *vct, Grid3DCU *grid){
#ifdef USEH5HUT

  H5input infile;

  infile.SetNameCycle(col->getinitfile(), col->getLast_cycle());

  infile.OpenFieldsFile("Node", nspec, col->getNxc()+1,
                                       col->getNyc()+1,
                                       col->getNzc()+1,
                                       vct->getCoordinates(),
                                       vct->getDims(),
                                       vct->getFieldComm());

  infile.ReadFields(EMf->getEx(), "Ex", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getEy(), "Ey", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getEz(), "Ez", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBx(), "Bx", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBy(), "By", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBz(), "Bz", grid->getNXN(), grid->getNYN(), grid->getNZN());

  for (int is = 0; is < nspec; is++){
    std::stringstream  ss;
    ss << is;
    std::string s_is = ss.str();
    infile.ReadFields(EMf->getRHOns(is), "rho_"+s_is, grid->getNXN(), grid->getNYN(), grid->getNZN());
  }

  infile.CloseFieldsFile();

  // initialize B on centers
    MPI_Barrier(MPIdata::get_PicGlobalComm());

  // Comm ghost nodes for B-field
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  grid->interpN2C(EMf->getBxc(), EMf->getBx());
  grid->interpN2C(EMf->getByc(), EMf->getBy());
  grid->interpN2C(EMf->getBzc(), EMf->getBz());

  // Comm ghost cells for B-field
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  // communicate E
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  for (int is = 0; is < nspec; is++)
    grid->interpN2C(EMf->getRHOcs(), is, EMf->getRHOns());

//---READ FROM THE CELLS:
//
//  infile.OpenFieldsFile("Cell", nspec, col->getNxc(),
//                                       col->getNyc(),
//                                       col->getNzc(),
//                                       vct->getCoordinates(),
//                                       vct->getDims(),
//                                       vct->getComm());
//
//  infile.ReadFields(EMf->getExc(), "Exc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getEyc(), "Eyc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getEzc(), "Ezc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getBxc(), "Bxc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getByc(), "Byc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getBzc(), "Bzc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//
//  for (int is = 0; is < nspec; is++){
//    std::stringstream  ss;
//    ss << is;
//    std::string s_is = ss.str();
//    infile.ReadFields(EMf->getRHOcs(is, 0), "rhoc_"+s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
//  }
//
//  infile.CloseFieldsFile();
//
//  // initialize B on nodes
//  grid->interpC2N(EMf->getBx(), EMf->getBxc());
//  grid->interpC2N(EMf->getBy(), EMf->getByc());
//  grid->interpC2N(EMf->getBz(), EMf->getBzc());
//
//  for (int is = 0; is < nspec; is++)
//    grid->interpC2N(EMf->getRHOns(), is, EMf->getRHOcs());
//
//---END READ FROM THE CELLS

#endif
}
#endif


template<typename T, int sz>
int size(T(&)[sz])
{
    return sz;
}

void WriteFieldsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & outputTag ,int cycle){

	//All VTK output at grid cells excluding ghost cells
	const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
	const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
	const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
	const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
	const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
	const int    nPoints = dimX*dimY*dimZ;
	MPI_File     fh;
	MPI_Status   status;

	if (outputTag.find("B", 0) != string::npos || outputTag.find("E", 0) != string::npos
			 || outputTag.find("Je", 0) != string::npos || outputTag.find("Ji", 0) != string::npos){

		const string tags0[]={"B", "E", "Je", "Ji"};
		float writebuffer[nzn-3][nyn-3][nxn-3][3];
		float tmpX, tmpY, tmpZ;

		for(int tagid=0;tagid<4;tagid++){
		 if (outputTag.find(tags0[tagid], 0) == string::npos) continue;

		 char   header[1024];
		 if (tags0[tagid].compare("B") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  tmpX = (float)EMf->getBxTot(ix+1, iy+1, iz+1);
						  tmpY = (float)EMf->getByTot(ix+1, iy+1, iz+1);
						  tmpZ = (float)EMf->getBzTot(ix+1, iy+1, iz+1);

						  writebuffer[iz][iy][ix][0] =  (fabs(tmpX) < 1E-16) ?0.0 : tmpX;
						  writebuffer[iz][iy][ix][1] =  (fabs(tmpY) < 1E-16) ?0.0 : tmpY;
						  writebuffer[iz][iy][ix][2] =  (fabs(tmpZ) < 1E-16) ?0.0 : tmpZ;
					  }

	      //Write VTK header
		  sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Magnetic Field from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d\n"
						   "VECTORS B float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (tags0[tagid].compare("E") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getEx(ix+1, iy+1, iz+1);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getEy(ix+1, iy+1, iz+1);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getEz(ix+1, iy+1, iz+1);
					  }

		  //Write VTK header
		   sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Electric Field from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS E float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (tags0[tagid].compare("Je") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 0);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 0);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 0);
					  }

		  //Write VTK header
		   sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Electron current from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS Je float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (tags0[tagid].compare("Ji") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 1);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 1);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 1);
					  }

		  //Write VTK header
		   sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Ion current from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS Ji float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);
		 }


		 if(EMf->isLittleEndian()){

			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][0],4);
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][1],4);
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][2],4);
					  }
		 }

		  int nelem = strlen(header);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_"<< tags0[tagid] << "_" << cycle << ".vtk";
		  MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

		  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
		  }

	      int err = MPI_File_set_view(fh, disp, EMf->getXYZeType(), EMf->getProcviewXYZ(), "native", MPI_INFO_NULL);
	      if(err){
	          dprintf("Error in MPI_File_set_view\n");

	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }

	      err = MPI_File_write_all(fh, writebuffer[0][0][0],3*(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
	      if(err){
		      int tcount=0;
		      MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			  dprintf(" wrote %i",  tcount);
	          dprintf("Error in write1\n");
	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }
	      MPI_File_close(&fh);
		}
	}

	if (outputTag.find("rho", 0) != string::npos){

		float writebufferRhoe[nzn-3][nyn-3][nxn-3];
		float writebufferRhoi[nzn-3][nyn-3][nxn-3];
		char   headerRhoe[1024];
		char   headerRhoi[1024];

		for(int iz=0;iz<nzn-3;iz++)
		  for(int iy=0;iy<nyn-3;iy++)
			  for(int ix= 0;ix<nxn-3;ix++){
				  writebufferRhoe[iz][iy][ix] = (float)EMf->getRHOns(ix+1, iy+1, iz+1, 0)*4*3.1415926535897;
				  writebufferRhoi[iz][iy][ix] = (float)EMf->getRHOns(ix+1, iy+1, iz+1, 1)*4*3.1415926535897;
			  }

		//Write VTK header
		sprintf(headerRhoe, "# vtk DataFile Version 2.0\n"
						   "Electron density from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "SCALARS rhoe float\n"
						   "LOOKUP_TABLE default\n",  dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		sprintf(headerRhoi, "# vtk DataFile Version 2.0\n"
						   "Ion density from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "SCALARS rhoi float\n"
						   "LOOKUP_TABLE default\n",  dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 if(EMf->isLittleEndian()){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  ByteSwap((unsigned char*) &writebufferRhoe[iz][iy][ix],4);
						  ByteSwap((unsigned char*) &writebufferRhoi[iz][iy][ix],4);
					  }
		 }

		  int nelem = strlen(headerRhoe);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_rhoe_" << cycle << ".vtk";
		  MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

		  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, headerRhoe, nelem, MPI_BYTE, &status);
		  }

	      int err = MPI_File_set_view(fh, disp, MPI_FLOAT, EMf->getProcview(), "native", MPI_INFO_NULL);
	      if(err){
	          dprintf("Error in MPI_File_set_view\n");

	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }

	      err = MPI_File_write_all(fh, writebufferRhoe[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
	      if(err){
		      int tcount=0;
		      MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			  dprintf(" wrote %i",  tcount);
	          dprintf("Error in write1\n");
	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }
	      MPI_File_close(&fh);

	      nelem = strlen(headerRhoi);
	      disp  = nelem*charsize;

	      filename.str("");
	      filename << col->getSaveDirName() << "/" << col->getSimName() << "_rhoi_" << cycle << ".vtk";
	      MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	      MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, headerRhoi, nelem, MPI_BYTE, &status);
		  }

	      err = MPI_File_set_view(fh, disp, MPI_FLOAT, EMf->getProcview(), "native", MPI_INFO_NULL);
		  if(err){
			  dprintf("Error in MPI_File_set_view\n");

			  int error_code = status.MPI_ERROR;
			  if (error_code != MPI_SUCCESS) {
				  char error_string[100];
				  int length_of_error_string, error_class;

				  MPI_Error_class(error_code, &error_class);
				  MPI_Error_string(error_class, error_string, &length_of_error_string);
				  dprintf("Error %s\n", error_string);
			  }
		  }

		  err = MPI_File_write_all(fh, writebufferRhoi[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
		  if(err){
			  int tcount=0;
			  MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			  dprintf(" wrote %i",  tcount);
			  dprintf("Error in write1\n");
			  int error_code = status.MPI_ERROR;
			  if (error_code != MPI_SUCCESS) {
				  char error_string[100];
				  int length_of_error_string, error_class;

				  MPI_Error_class(error_code, &error_class);
				  MPI_Error_string(error_class, error_string, &length_of_error_string);
				  dprintf("Error %s\n", error_string);
			  }
		  }
		  MPI_File_close(&fh);
	}
}

void WriteMomentsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & outputTag ,int cycle, float*** momentswritebuffer){

	//All VTK output at grid cells excluding ghost cells
	const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
	const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
	const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
	const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
	const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
	const int    nPoints = dimX*dimY*dimZ;
	MPI_File     fh;
	MPI_Status   status;
	const string momentstags[]={"rho", "PXX", "PXY", "PXZ", "PYY", "PYZ", "PZZ"};
	const int    tagsize = size(momentstags);
	const string outputtag = col->getMomentsOutputTag();
	const int ns = col->getNs();

	for(int tagid=0; tagid<tagsize; tagid++){
		 if (outputtag.find(momentstags[tagid], 0) == string::npos) continue;

		 for(int si=0;si<ns;si++){
			 char  header[1024];
			 if (momentstags[tagid].compare("rho") == 0){
				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getRHOns(ix+1, iy+1, iz+1, si)*4*3.1415926535897;

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s%d density from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
					"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"rhoe":"rhoi");
		 }else if(momentstags[tagid].compare("PXX") == 0){
				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
				    for(int ix= 0;ix<nxn-3;ix++)
				      momentswritebuffer[iz][iy][ix] = (float)EMf->getpXXsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
						"%s%d pressure PXX from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
					"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PXXe":"PXXi");
		}else if(momentstags[tagid].compare("PXY") == 0){

			for(int iz=0;iz<nzn-3;iz++)
			  for(int iy=0;iy<nyn-3;iy++)
				  for(int ix= 0;ix<nxn-3;ix++)
					  momentswritebuffer[iz][iy][ix] = (float)EMf->getpXYsn(ix+1, iy+1, iz+1, si);

			//Write VTK header
			sprintf(header, "# vtk DataFile Version 2.0\n"
					"%s%d pressure PXY from iPIC3D\n"
							   "BINARY\n"
							   "DATASET STRUCTURED_POINTS\n"
							   "DIMENSIONS %d %d %d\n"
							   "ORIGIN 0 0 0\n"
							   "SPACING %f %f %f\n"
							   "POINT_DATA %d \n"
							   "SCALARS %s float\n"
				"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PXYe":"PXYi");
		}else if(momentstags[tagid].compare("PXZ") == 0){

			for(int iz=0;iz<nzn-3;iz++)
			  for(int iy=0;iy<nyn-3;iy++)
				  for(int ix= 0;ix<nxn-3;ix++)
					  momentswritebuffer[iz][iy][ix] = (float)EMf->getpXZsn(ix+1, iy+1, iz+1, si);

			//Write VTK header
			sprintf(header, "# vtk DataFile Version 2.0\n"
							   "%s%d pressure PXZ from iPIC3D\n"
							   "BINARY\n"
							   "DATASET STRUCTURED_POINTS\n"
							   "DIMENSIONS %d %d %d\n"
							   "ORIGIN 0 0 0\n"
							   "SPACING %f %f %f\n"
							   "POINT_DATA %d \n"
							   "SCALARS %s float\n"
				"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PXZe":"PXZi");
		}else if(momentstags[tagid].compare("PYY") == 0){

			for(int iz=0;iz<nzn-3;iz++)
			  for(int iy=0;iy<nyn-3;iy++)
				  for(int ix= 0;ix<nxn-3;ix++)
					  momentswritebuffer[iz][iy][ix] = (float)EMf->getpYYsn(ix+1, iy+1, iz+1, si);

			//Write VTK header
			sprintf(header, "# vtk DataFile Version 2.0\n"
							   "%s%d pressure PYY from iPIC3D\n"
							   "BINARY\n"
							   "DATASET STRUCTURED_POINTS\n"
							   "DIMENSIONS %d %d %d\n"
							   "ORIGIN 0 0 0\n"
							   "SPACING %f %f %f\n"
							   "POINT_DATA %d \n"
							   "SCALARS %s float\n"
				"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PYYe":"PYYi");
		}else if(momentstags[tagid].compare("PYZ") == 0){

			for(int iz=0;iz<nzn-3;iz++)
			  for(int iy=0;iy<nyn-3;iy++)
				  for(int ix= 0;ix<nxn-3;ix++)
					  momentswritebuffer[iz][iy][ix] = (float)EMf->getpYZsn(ix+1, iy+1, iz+1, si);

			//Write VTK header
			sprintf(header, "# vtk DataFile Version 2.0\n"
							   "%s%d pressure PYZ from iPIC3D\n"
							   "BINARY\n"
							   "DATASET STRUCTURED_POINTS\n"
							   "DIMENSIONS %d %d %d\n"
							   "ORIGIN 0 0 0\n"
							   "SPACING %f %f %f\n"
							   "POINT_DATA %d \n"
							   "SCALARS %s float\n"
				"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PYZe":"PYZi");
		}else if(momentstags[tagid].compare("PZZ") == 0){

			for(int iz=0;iz<nzn-3;iz++)
			  for(int iy=0;iy<nyn-3;iy++)
				  for(int ix= 0;ix<nxn-3;ix++)
					  momentswritebuffer[iz][iy][ix] = (float)EMf->getpZZsn(ix+1, iy+1, iz+1, si);

			//Write VTK header
			sprintf(header, "# vtk DataFile Version 2.0\n"
							   "%s%d pressure PZZ from iPIC3D\n"
							   "BINARY\n"
							   "DATASET STRUCTURED_POINTS\n"
							   "DIMENSIONS %d %d %d\n"
							   "ORIGIN 0 0 0\n"
							   "SPACING %f %f %f\n"
							   "POINT_DATA %d \n"
							   "SCALARS %s float\n"
				"LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PZZe":"PZZi");
		}

		 if(EMf->isLittleEndian()){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  ByteSwap((unsigned char*) &momentswritebuffer[iz][iy][ix],4);
					  }
		 }

		  int nelem = strlen(header);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_" << momentstags[tagid] << ((si%2==0)?"e":"i")<< si  << "_" << cycle << ".vtk";
		  MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
		  }

	      int error_code = MPI_File_set_view(fh, disp, MPI_FLOAT, EMf->getProcview(), "native", MPI_INFO_NULL);
	      if (error_code != MPI_SUCCESS) {
			char error_string[100];
			int length_of_error_string, error_class;
			MPI_Error_class(error_code, &error_class);
			MPI_Error_string(error_class, error_string, &length_of_error_string);
			dprintf("Error in MPI_File_set_view: %s\n", error_string);
		  }

	      error_code = MPI_File_write_all(fh, momentswritebuffer[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
	      if(error_code != MPI_SUCCESS){
		      int tcount=0;
		      MPI_Get_count(&status, MPI_FLOAT, &tcount);
	          char error_string[100];
	          int length_of_error_string, error_class;
	          MPI_Error_class(error_code, &error_class);
	          MPI_Error_string(error_class, error_string, &length_of_error_string);
	          dprintf("Error in MPI_File_write_all: %s, wrote %d MPI_FLOAT\n", error_string,tcount);
	      }
	      MPI_File_close(&fh);
		 }//END OF SPECIES
	}//END OF TAGS
}

void WriteFieldsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & outputTag ,int cycle,float**** fieldwritebuffer){

  //All VTK output at grid cells excluding ghost cells
  const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
  const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
  const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
  const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
  const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
  const int    nPoints = dimX*dimY*dimZ;
  MPI_File     fh;
  MPI_Status   status;
  const string fieldtags[]={"B", "E", "Je", "Ji","Je2", "Ji3"};
  const int    tagsize = size(fieldtags);
  const string outputtag = col->getFieldOutputTag();

  for(int tagid=0; tagid<tagsize; tagid++){

   if (outputTag.find(fieldtags[tagid], 0) == string::npos) continue;

   char   header[1024];
   if (fieldtags[tagid].compare("B") == 0){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            fieldwritebuffer[iz][iy][ix][0] =  (float)EMf->getBxTot(ix+1, iy+1, iz+1);
            fieldwritebuffer[iz][iy][ix][1] =  (float)EMf->getByTot(ix+1, iy+1, iz+1);
            fieldwritebuffer[iz][iy][ix][2] =  (float)EMf->getBzTot(ix+1, iy+1, iz+1);
          }

     //Write VTK header
     sprintf(header, "# vtk DataFile Version 2.0\n"
             "Magnetic Field from iPIC3D\n"
             "BINARY\n"
             "DATASET STRUCTURED_POINTS\n"
             "DIMENSIONS %d %d %d\n"
             "ORIGIN 0 0 0\n"
             "SPACING %f %f %f\n"
             "POINT_DATA %d\n"
             "VECTORS B float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

   }else if (fieldtags[tagid].compare("E") == 0){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getEx(ix+1, iy+1, iz+1);
            fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getEy(ix+1, iy+1, iz+1);
            fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getEz(ix+1, iy+1, iz+1);
          }

     //Write VTK header
     sprintf(header, "# vtk DataFile Version 2.0\n"
             "Electric Field from iPIC3D\n"
             "BINARY\n"
             "DATASET STRUCTURED_POINTS\n"
             "DIMENSIONS %d %d %d\n"
             "ORIGIN 0 0 0\n"
             "SPACING %f %f %f\n"
             "POINT_DATA %d \n"
             "VECTORS E float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

   }else if (fieldtags[tagid].compare("Je") == 0){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 0);
            fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 0);
            fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 0);
          }

     //Write VTK header
     sprintf(header, "# vtk DataFile Version 2.0\n"
             "Electron current from iPIC3D\n"
             "BINARY\n"
             "DATASET STRUCTURED_POINTS\n"
             "DIMENSIONS %d %d %d\n"
             "ORIGIN 0 0 0\n"
             "SPACING %f %f %f\n"
             "POINT_DATA %d \n"
             "VECTORS Je float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

   }else if (fieldtags[tagid].compare("Ji") == 0){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 1);
            fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 1);
            fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 1);
          }

     //Write VTK header
     sprintf(header, "# vtk DataFile Version 2.0\n"
             "Ion current from iPIC3D\n"
             "BINARY\n"
             "DATASET STRUCTURED_POINTS\n"
             "DIMENSIONS %d %d %d\n"
             "ORIGIN 0 0 0\n"
             "SPACING %f %f %f\n"
             "POINT_DATA %d \n"
             "VECTORS Ji float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);
   }else if (fieldtags[tagid].compare("Je2") == 0){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 2);
            fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 2);
            fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 2);
          }

     //Write VTK header
     sprintf(header, "# vtk DataFile Version 2.0\n"
         "Electron2 current from iPIC3D\n"
             "BINARY\n"
             "DATASET STRUCTURED_POINTS\n"
             "DIMENSIONS %d %d %d\n"
             "ORIGIN 0 0 0\n"
             "SPACING %f %f %f\n"
             "POINT_DATA %d \n"
             "VECTORS Je float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

   }else if (fieldtags[tagid].compare("Ji3") == 0){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 3);
            fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 3);
            fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 3);
          }

     //Write VTK header
     sprintf(header, "# vtk DataFile Version 2.0\n"
         "Ion3 current from iPIC3D\n"
             "BINARY\n"
             "DATASET STRUCTURED_POINTS\n"
             "DIMENSIONS %d %d %d\n"
             "ORIGIN 0 0 0\n"
             "SPACING %f %f %f\n"
             "POINT_DATA %d \n"
             "VECTORS Ji float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);
   }

   if(EMf->isLittleEndian()){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
            ByteSwap((unsigned char*) &fieldwritebuffer[iz][iy][ix][0],4);
            ByteSwap((unsigned char*) &fieldwritebuffer[iz][iy][ix][1],4);
            ByteSwap((unsigned char*) &fieldwritebuffer[iz][iy][ix][2],4);
          }
   }

    int nelem = strlen(header);
    int charsize=sizeof(char);
    MPI_Offset disp = nelem*charsize;

    ostringstream filename;
    filename << col->getSaveDirName() << "/" << col->getSimName() << "_"<< fieldtags[tagid] << "_" << cycle << ".vtk";
    MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    if (vct->getCartesian_rank()==0){
      MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
    }

      int error_code = MPI_File_set_view(fh, disp, EMf->getXYZeType(), EMf->getProcviewXYZ(), "native", MPI_INFO_NULL);
      if (error_code != MPI_SUCCESS) {
    char error_string[100];
    int length_of_error_string, error_class;

    MPI_Error_class(error_code, &error_class);
    MPI_Error_string(error_class, error_string, &length_of_error_string);
    dprintf("Error in MPI_File_set_view: %s\n", error_string);
    }

      error_code = MPI_File_write_all(fh, fieldwritebuffer[0][0][0],(nxn-3)*(nyn-3)*(nzn-3),EMf->getXYZeType(), &status);
      if(error_code != MPI_SUCCESS){
        int tcount=0;
        MPI_Get_count(&status, EMf->getXYZeType(), &tcount);
          char error_string[100];
          int length_of_error_string, error_class;
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          dprintf("Error in MPI_File_write_all: %s, wrote %d EMf->getXYZeType()\n", error_string,tcount);
      }
      MPI_File_close(&fh);
  }
}

void WriteTemperatureVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & outputTag ,int cycle,float**** temperaturewritebuffer){

  //All VTK output at grid cells excluding ghost cells
  const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
  const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
  const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
  const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
  const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
  const int    nPoints = dimX*dimY*dimZ;
  MPI_File     fh;
  MPI_Status   status;
  const string temperaturetags[]={"Tcart","Tperpar"};
  const int    tagsize = size(temperaturetags);
  const string outputtag = col->getTemperatureOutputTag();
  const int ns = col->getNs();

  double Tcar[nzn-3][nyn-3][nxn-3][6] = {};
  double Tperpar[nzn-3][nyn-3][nxn-3][6] = {};

  for(int si=0;si<ns;si++){
   char   header[1024];

     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++){
              double pxx = (float)EMf->getpXXsn(ix+1, iy+1, iz+1, si);
              double pxy = (float)EMf->getpXYsn(ix+1, iy+1, iz+1, si);
              double pxz = (float)EMf->getpXZsn(ix+1, iy+1, iz+1, si);
              double pyy = (float)EMf->getpYYsn(ix+1, iy+1, iz+1, si);
              double pyz = (float)EMf->getpYZsn(ix+1, iy+1, iz+1, si);
              double pzz = (float)EMf->getpZZsn(ix+1, iy+1, iz+1, si);

              double rho = (float)EMf->getRHOns(ix+1, iy+1, iz+1, si);
              if (rho==0){
                rho=0.00001;
              }

              double Jx = (float)EMf->getJxs(ix+1, iy+1, iz+1, si);
              double Jy = (float)EMf->getJys(ix+1, iy+1, iz+1, si);
              double Jz = (float)EMf->getJzs(ix+1, iy+1, iz+1, si);

              double qom = col->getQOM(si);

              double Bx = (float)EMf->getBxTot(ix+1, iy+1, iz+1);
              double By = (float)EMf->getByTot(ix+1, iy+1, iz+1);
              double Bz = (float)EMf->getBzTot(ix+1, iy+1, iz+1);
              double modB = sqrt(pow(Bx,2) + pow(By,2) + pow(Bz,2));

              double a1 = Bx/modB;
              double a2 = By/modB;
              double a3 = Bz/modB;

              double b1;
              double b2;
              double b3;
              if (a2 != 0 || a3 != 0){
                b1 = 0;
                b2 = a3;
                b3 = -a2;
              }else{
                b1 = -a3;
                b2 = 0;
                b3 = a1;
              }
              double modb = sqrt(pow(b1,2) + pow(b2,2) + pow(b3,2));
              b1 = b1/modb;
              b2 = b2/modb;
              b3 = b3/modb;

              double c1 = a2*b3 - a3*b2;
              double c2 = a3*b1 - a1*b3;
              double c3 = a1*b2 - a2*b1;
              
              double Txx = (pxx-((Jx*Jx)/rho))/(rho*fabs(qom));
              double Txy = (pxy-((Jx*Jy)/rho))/(rho*fabs(qom));
              double Txz = (pxz-((Jx*Jz)/rho))/(rho*fabs(qom));
              double Tyy = (pyy-((Jy*Jy)/rho))/(rho*fabs(qom));
              double Tyz = (pyz-((Jy*Jz)/rho))/(rho*fabs(qom));
              double Tzz = (pzz-((Jz*Jz)/rho))/(rho*fabs(qom));

              Tcar[iz][iy][ix][0] = Txx;
              Tcar[iz][iy][ix][1] = Txy;
              Tcar[iz][iy][ix][2] = Txz;
              Tcar[iz][iy][ix][3] = Tyy;
              Tcar[iz][iy][ix][4] = Tyz;
              Tcar[iz][iy][ix][5] = Tzz;

              double Tai = a1*b1*Txx + a2*b2*Tyy + a3*b3*Tzz + (a1*b2 + a2*b1)*Txy + (a1*b3 + a3*b1)*Txz + (a2*b3 + a3*b2)*Tyz;
              double Tbi = a1*c1*Txx + a2*c2*Tyy + a3*c3*Tzz + (a1*c2 + a2*c1)*Txy + (a1*c3 + a3*c1)*Txz + (a2*c3 + a3*c2)*Tyz;
              double Tci = b1*c1*Txx + b2*c2*Tyy + b3*c3*Tzz + (b1*c2 + b2*c1)*Txy + (b1*c3 + b3*c1)*Txz + (b2*c3 + b3*c2)*Tyz;
              double Tpar = a1*a1*Txx + a2*a2*Tyy + a3*a3*Tzz + 2*(a1*a2*Txy + a1*a3*Txz + a2*a3*Tyz);
              double Tp1i = b1*b1*Txx + b2*b2*Tyy + b3*b3*Tzz + 2*(b1*b2*Txy + b1*b3*Txz + b2*b3*Tyz);
              double Tp2i = c1*c1*Txx + c2*c2*Tyy + c3*c3*Tzz + 2*(c1*c2*Txy + c1*c3*Txz + c2*c3*Tyz);

              double theta = 0.5*atan((Tp2i-Tp1i)/(2*Tci));

              double Ta = cos(theta)*Tai + sin(theta)*Tbi;
              double Tb = cos(theta)*Tbi - sin(theta)*Tai;
              double Tc = cos(theta)*cos(theta)*Tci - sin(theta)*sin(theta)*Tci + cos(theta)*sin(theta)*(Tp2i-Tp1i); //cos(2*theta)*Tci + 0.5*sin(2*theta)*(Tp2i - Tp1i);
              double Tperp1 = 2*cos(theta)*sin(theta)*Tci + cos(theta)*cos(theta)*Tp1i + sin(theta)*sin(theta)*Tp2i; //Tci + 0.5*sin(2*theta)*(Tp2i + Tp1i);
              double Tperp2 = -2*cos(theta)*sin(theta)*Tci + cos(theta)*cos(theta)*Tp2i + sin(theta)*sin(theta)*Tp1i;

              Tperpar[iz][iy][ix][0] = Tpar;
              Tperpar[iz][iy][ix][1] = Tperp1;
              Tperpar[iz][iy][ix][2] = Tperp2;
              Tperpar[iz][iy][ix][3] = Ta;
              Tperpar[iz][iy][ix][4] = Tb;
              Tperpar[iz][iy][ix][5] = Tc;
          }

   for(int tagid=0; tagid<tagsize; tagid++){

   if (outputTag.find(temperaturetags[tagid], 0) == string::npos) continue;

   if (temperaturetags[tagid].compare("Tcart") == 0){
    //Write VTK header
      sprintf(header, "# vtk DataFile Version 2.0\n"
          "%s%d temperature tensor in cartesian coordinates from iPIC3D\n"
                 "BINARY\n"
                 "DATASET STRUCTURED_POINTS\n"
                 "DIMENSIONS %d %d %d\n"
                 "ORIGIN 0 0 0\n"
                 "SPACING %f %f %f\n"
                 "POINT_DATA %d \n"
                 "TENSORS6 %s float\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"Tcart_e":"Tcart_i");

      for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++)
              for(int s=0;s<6;s++){
                temperaturewritebuffer[iz][iy][ix][s] = Tcar[iz][iy][ix][s];
              }

   }else if(temperaturetags[tagid].compare("Tperpar") == 0){
    //Write VTK header
      sprintf(header, "# vtk DataFile Version 2.0\n"
          "%s%d temperature tensor in coordinates relative to the magnetic field from iPIC3D\n"
                 "BINARY\n"
                 "DATASET STRUCTURED_POINTS\n"
                 "DIMENSIONS %d %d %d\n"
                 "ORIGIN 0 0 0\n"
                 "SPACING %f %f %f\n"
                 "POINT_DATA %d \n"
                 "TENSORS6 %s float\n",(si%2==0)?"Electron":"Ion",si,dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"Tperpar_e":"Tperpar_i");

      for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++)
              for(int s=0;s<6;s++){
                temperaturewritebuffer[iz][iy][ix][s] = Tperpar[iz][iy][ix][s];
              }

   }

   if(EMf->isLittleEndian()){
     for(int iz=0;iz<nzn-3;iz++)
        for(int iy=0;iy<nyn-3;iy++)
          for(int ix= 0;ix<nxn-3;ix++)
            for(int s=0;s<6;s++){
              ByteSwap((unsigned char*) &temperaturewritebuffer[iz][iy][ix][s],4);
          }
   }

    int nelem = strlen(header);
    int charsize=sizeof(char);
    MPI_Offset disp = nelem*charsize;


    ostringstream filename;
    filename << col->getSaveDirName() << "/" << col->getSimName() << "_"<< temperaturetags[tagid] << ((si%2==0)?"e":"i")<< si  << "_" << cycle << ".vtk";
    MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    if (vct->getCartesian_rank()==0){
      MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
    }

    int error_code = MPI_File_set_view(fh, disp, EMf->getTensorType(), EMf->getProcviewTensor(), "native", MPI_INFO_NULL);
    if (error_code != MPI_SUCCESS) {
      char error_string[100];
      int length_of_error_string, error_class;

      MPI_Error_class(error_code, &error_class);
      MPI_Error_string(error_class, error_string, &length_of_error_string);
      dprintf("Error in MPI_File_set_view: %s\n", error_string);
    }

      error_code = MPI_File_write_all(fh, temperaturewritebuffer[0][0][0],(nxn-3)*(nyn-3)*(nzn-3),EMf->getTensorType(), &status);
      if(error_code != MPI_SUCCESS){
        int tcount=0;
        MPI_Get_count(&status, EMf->getTensorType(), &tcount);
          char error_string[100];
          int length_of_error_string, error_class;
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          dprintf("Error in MPI_File_write_all: %s, wrote %d EMf->getTensorType()\n", error_string,tcount);
      }
      MPI_File_close(&fh);
    }//END OF TAGS
  }//END OF SPECIES
}

void WriteSpectraVTK(Grid3DCU *grid, Particles3D *part, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & outputTag ,int cycle, float**** spectrawritebuffere, float**** spectrawritebufferi){

  //All VTK output at grid cells excluding ghost cells
  const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
  const int nxc  =grid->getNXC(),nyc  = grid->getNYC(),nzc =grid->getNZC();
  const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
  const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
  const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
  const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
  const double dx = col->getDx(), dy = col->getDy(), dz = col->getDz();
  const int    DeltaX = col->getDeltaX(), DeltaY = col->getDeltaY(), DeltaZ = col->getDeltaZ();
  const double Estarti = col->getEstarti(), Eendi = col->getEendi(), dEi = col->getdEi();
  const double Estarte = col->getEstarte(), Eende = col->getEende(), dEe = col->getdEe();
  const double xWidth = dx*(nxn-3), yWidth = dy*(nyn-3), zWidth = dz*(nzn-3);
  const double xStart = vct->getCoordinates(0) * xWidth, yStart = vct->getCoordinates(1) * yWidth, zStart = vct->getCoordinates(2) * zWidth;
  const int    nPoints = dimX*dimY*dimZ/(DeltaX*DeltaY*DeltaZ);
  MPI_File     fh;
  MPI_Status   status;
  const string spectratags[]={"Stot", "Spar", "Sperp"};
  const int    tagsize = size(spectratags);
  const string outputtag = col->getSpectraOutputTag();
  const int ns = col->getNs();
  const int sizeX = (nxn-3)/DeltaX, sizeY = (nyn-3)/DeltaY, sizeZ = (nzn-3)/DeltaZ;
  double Estart = 0, Eend = 0, dE = 0;

  for(int si=0;si<ns;si++){

      //si%2==0 for electrons, si%2==1 for ions.
      if (si%2==0){
        Estart = Estarte;
        Eend = Eende;
        dE = dEe;
      }else{
        Estart = Estarti;
        Eend = Eendi;
        dE = dEi;
      }

      int Nspec = (Eend-Estart)/dE+1.001;
      double spectra[Nspec];
      for(int s=0;s<Nspec;s++){
        spectra[s] = pow(10,Estart + dE*s);
      }

      if (vct->getCartesian_rank()==0){
        FILE *fptr;
        ostringstream savedir;
        savedir << col->getSaveDirName() << "/EnergyBins_" << ((si%2==0)?"Electron":"Ion") << ".txt"; 
        fptr = fopen(savedir.str().c_str(),"w");

        for(int s=0;s<Nspec;s++){
          fprintf(fptr,"%e\n",spectra[s]);
        }
        fclose(fptr);
      }

      double Stot[sizeZ][sizeY][sizeX][Nspec+1] = {};
      double Spar[sizeZ][sizeY][sizeX][Nspec+1] = {};
      double Sperp[sizeZ][sizeY][sizeX][Nspec+1] = {};

       double qom = col->getQOM(si);

       if (si%2==0){
       for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebuffere[iz][iy][ix][s] = 0;
                Stot[iz][iy][ix][s] = 0;
                Spar[iz][iy][ix][s] = 0;
                Sperp[iz][iy][ix][s] = 0;
              }
        }else{
          for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebufferi[iz][iy][ix][s] = 0;
                Stot[iz][iy][ix][s] = 0;
                Spar[iz][iy][ix][s] = 0;
                Sperp[iz][iy][ix][s] = 0;
              }
        }

       char  header[1024];

        for (int n = 0; n < part[si].getNOP(); n++){
          double x = part[si].getX(n);
          double y = part[si].getY(n);
          double z = part[si].getZ(n);
          double u = part[si].getU(n);
          double v = part[si].getV(n);
          double w = part[si].getW(n);

          for(int iz=0;iz<sizeZ;iz++){
            if(z > zStart + dz*DeltaZ*(iz+1)) continue;

            for(int iy=0;iy<sizeY;iy++){
              if(y > yStart + dy*DeltaY*(iy+1)) continue;

              for(int ix=0;ix<sizeX;ix++){
                if(x > xStart + dx*DeltaX*(ix+1)) continue;

                double Bx = 0;
                double By = 0;
                double Bz = 0;

                for(int cz=0;cz<DeltaZ;cz++)
                  for(int cy=0;cy<DeltaY;cy++)
                    for(int cx=0;cx<DeltaX;cx++){
                      double Bxi = (float)EMf->getBxTot(ix+cx+1, iy+cy+1, iz+cz+1);
                      double Byi = (float)EMf->getByTot(ix+cx+1, iy+cy+1, iz+cz+1);
                      double Bzi = (float)EMf->getBzTot(ix+cx+1, iy+cy+1, iz+cz+1);
                      Bx += Bxi;
                      By += Byi;
                      Bz += Bzi;
                    }
                double modB = sqrt(pow(Bx,2) + pow(By,2) + pow(Bz,2));

                ///////// FIRST COORDINATE CHANGE
                double a1 = Bx/modB;
                double a2 = By/modB;
                double a3 = Bz/modB;

                double b1;
                double b2;
                double b3;
                if (a2 != 0 || a3 != 0){
                  b1 = 0;
                  b2 = a3;
                  b3 = -a2;
                }else{
                  b1 = -a3;
                  b2 = 0;
                  b3 = a1;
                }
                double modb = sqrt(pow(b1,2) + pow(b2,2) + pow(b3,2));
                b1 = b1/modb;
                b2 = b2/modb;
                b3 = b3/modb;

                double c1 = a2*b3 - a3*b2;
                double c2 = a3*b1 - a1*b3;
                double c3 = a1*b2 - a2*b1;

                double vpar = a1*u + a2*v + a3*w;
                double vperp1 = b1*u + b2*v + b3*w;
                double vperp2 = c1*u + c2*v + c3*w;

                ///////// SPECTRA CALC
                double Epar = 0.5*vpar*vpar/fabs(qom);
                double Eperp1 = 0.5*vperp1*vperp1/fabs(qom);
                double Eperp2 = 0.5*vperp2*vperp2/fabs(qom);

                for(int s=0;s<Nspec+1;s++){
                  if(s==Nspec){
                    if(Epar+Eperp1+Eperp2 > spectra[s-1]){
                      Stot[iz][iy][ix][s] += fabs(part[si].getQ(n));
                    }
                  }else if(Epar+Eperp1+Eperp2 < spectra[s]){
                    Stot[iz][iy][ix][s] += fabs(part[si].getQ(n));
                    break;
                  }else{
                    continue;
                  }
                }
                for(int s=0;s<Nspec+1;s++){
                  if(s==Nspec){
                    if(Epar > spectra[s-1]){
                      Spar[iz][iy][ix][s] += fabs(part[si].getQ(n));
                    }
                  }else if(Epar < spectra[s]){
                    Spar[iz][iy][ix][s] += fabs(part[si].getQ(n));
                    break;
                  }else{
                    continue;
                  }
                }
                for(int s=0;s<Nspec+1;s++){
                  if(s==Nspec){
                    if(Eperp1+Eperp2 > spectra[s-1]){
                      Sperp[iz][iy][ix][s] += fabs(part[si].getQ(n));
                    }
                  }else if(Eperp1+Eperp2 < spectra[s]){
                    Sperp[iz][iy][ix][s] += fabs(part[si].getQ(n));
                    break;
                  }else{
                    continue;
                  }
                }
                break;
              }
              break;
            }
            break;
          }
        }

    for(int tagid=0; tagid<tagsize; tagid++){
      if (outputtag.find(spectratags[tagid], 0) == string::npos) continue;

     
      if(spectratags[tagid].compare("Stot") == 0){
        //Write VTK header
        sprintf(header, "# vtk DataFile Version 2.0\n"
            "%s%d total energy spectra from iPIC3D\n"
                   "BINARY\n"
                   "DATASET STRUCTURED_POINTS\n"
                   "DIMENSIONS %d %d %d\n"
                   "ORIGIN 0 0 0\n"
                   "SPACING %f %f %f\n"
                   "POINT_DATA %d \n"
                   "SCALARS %s float %d \n"
          "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX/DeltaX,dimY/DeltaY,dimZ/DeltaZ, spaceX*DeltaX,spaceY*DeltaY,spaceZ*DeltaZ, nPoints,(si%2==0)?"Se":"Si", Nspec+1);

      if (si%2==0){
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebuffere[iz][iy][ix][s] = Stot[iz][iy][ix][s];
              }
      }else{
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebufferi[iz][iy][ix][s] = Stot[iz][iy][ix][s];
              }
      }
    }

    if(spectratags[tagid].compare("Spar") == 0){
        //Write VTK header
        sprintf(header, "# vtk DataFile Version 2.0\n"
            "%s%d parallel energy spectra from iPIC3D\n"
                   "BINARY\n"
                   "DATASET STRUCTURED_POINTS\n"
                   "DIMENSIONS %d %d %d\n"
                   "ORIGIN 0 0 0\n"
                   "SPACING %f %f %f\n"
                   "POINT_DATA %d \n"
                   "SCALARS %s float %d \n"
          "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX/DeltaX,dimY/DeltaY,dimZ/DeltaZ, spaceX*DeltaX,spaceY*DeltaY,spaceZ*DeltaZ, nPoints,(si%2==0)?"Se":"Si", Nspec+1);


      if (si%2==0){
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebuffere[iz][iy][ix][s] = Spar[iz][iy][ix][s];
              }
      }else{
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebufferi[iz][iy][ix][s] = Spar[iz][iy][ix][s];
              }
      }
    }

    if(spectratags[tagid].compare("Sperp") == 0){
        //Write VTK header
        sprintf(header, "# vtk DataFile Version 2.0\n"
            "%s%d perpendicular energy spectra from iPIC3D\n"
                   "BINARY\n"
                   "DATASET STRUCTURED_POINTS\n"
                   "DIMENSIONS %d %d %d\n"
                   "ORIGIN 0 0 0\n"
                   "SPACING %f %f %f\n"
                   "POINT_DATA %d \n"
                   "SCALARS %s float %d \n"
          "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",si,dimX/DeltaX,dimY/DeltaY,dimZ/DeltaZ, spaceX*DeltaX,spaceY*DeltaY,spaceZ*DeltaZ, nPoints,(si%2==0)?"Se":"Si", Nspec+1);


      if (si%2==0){
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebuffere[iz][iy][ix][s] = Sperp[iz][iy][ix][s];
              }
      }else{
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++)
              for(int s=0;s<Nspec+1;s++){
                spectrawritebufferi[iz][iy][ix][s] = Sperp[iz][iy][ix][s];
              }
      }
    }


     if(EMf->isLittleEndian()){
      if (si%2==0){
       for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++){
              for(int s=0;s<Nspec+1;s++){
                ByteSwap((unsigned char*) &spectrawritebuffere[iz][iy][ix][s],4);
              }
            }
      }else{
        for(int iz=0;iz<sizeZ;iz++)
          for(int iy=0;iy<sizeY;iy++)
            for(int ix= 0;ix<sizeX;ix++){
              for(int s=0;s<Nspec+1;s++){
                ByteSwap((unsigned char*) &spectrawritebufferi[iz][iy][ix][s],4);
              }
            }
      }
     }

      int nelem = strlen(header);
      int charsize=sizeof(char);
      MPI_Offset disp = nelem*charsize;

      ostringstream filename;
      filename << col->getSaveDirName() << "/" << col->getSimName() << "_" << spectratags[tagid] << ((si%2==0)?"e":"i")<< si  << "_" << cycle << ".vtk";
      MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

      if (vct->getCartesian_rank()==0){
        MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
      }

      MPI_Datatype  spectracomp, spectraview;

       //create process file view
      int  size[3], subsize[3], start[3];

      //3D subarray - reverse X, Z
      subsize[0] = sizeZ; subsize[1] = sizeY; subsize[2] = sizeX;
      size[0] = sizeZ*vct->getZLEN();size[1] = sizeY*vct->getYLEN();size[2] = sizeX*vct->getXLEN();
      start[0]= vct->getCoordinates(2)*subsize[0];
      start[1]= vct->getCoordinates(1)*subsize[1];
      start[2]= vct->getCoordinates(0)*subsize[2];

      MPI_Type_contiguous(Nspec+1,MPI_FLOAT, &spectracomp);
      MPI_Type_commit(&spectracomp);

      MPI_Type_create_subarray(3, size, subsize, start,MPI_ORDER_C, spectracomp, &spectraview);
      MPI_Type_commit(&spectraview);

      int error_code = MPI_File_set_view(fh, disp, spectracomp, spectraview, "native", MPI_INFO_NULL);
      if (error_code != MPI_SUCCESS) {
        char error_string[100];
        int length_of_error_string, error_class;

        MPI_Error_class(error_code, &error_class);
        MPI_Error_string(error_class, error_string, &length_of_error_string);
        dprintf("Error in MPI_File_set_view: %s\n", error_string);
      }
        if (si%2==0){
          error_code = MPI_File_write_all(fh, spectrawritebuffere[0][0][0],sizeX*sizeY*sizeZ,spectracomp, &status);
        }else{
          error_code = MPI_File_write_all(fh, spectrawritebufferi[0][0][0],sizeX*sizeY*sizeZ,spectracomp, &status);
        }
        if(error_code != MPI_SUCCESS){
          int tcount=0;
          MPI_Get_count(&status, spectracomp, &tcount);
            char error_string[100];
            int length_of_error_string, error_class;
            MPI_Error_class(error_code, &error_class);
            MPI_Error_string(error_class, error_string, &length_of_error_string);
            dprintf("Error in MPI_File_write_all: %s, wrote %d spectracomp\n", error_string,tcount);
        }
        MPI_File_close(&fh);
        MPI_Type_free(&spectraview);
        MPI_Type_free(&spectracomp);
     }//END OF TAGS
  }//END OF SPECIES
}

int WriteFieldsVTKNonblk(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct,int cycle,
		float**** fieldwritebuffer,MPI_Request requestArr[4],MPI_File fhArr[4]){

	//All VTK output at grid cells excluding ghost cells
	const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
	const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
	const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
	const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
	const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
	const int    nPoints = dimX*dimY*dimZ;
	const string fieldtags[]={"B", "E", "Je", "Ji"};
	const int fieldtagsize = size(fieldtags);
	const string fieldoutputtag = col->getFieldOutputTag();
	int counter=0,error_code;

	for(int tagid=0; tagid<fieldtagsize; tagid++){
		 if (fieldoutputtag.find(fieldtags[tagid], 0) == string::npos) continue;

		 fieldwritebuffer = &(fieldwritebuffer[counter]);
		 char  header[1024];
		 if (fieldtags[tagid].compare("B") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  fieldwritebuffer[iz][iy][ix][0] =  (float)EMf->getBxTot(ix+1, iy+1, iz+1);
						  fieldwritebuffer[iz][iy][ix][1] =  (float)EMf->getByTot(ix+1, iy+1, iz+1);
						  fieldwritebuffer[iz][iy][ix][2] =  (float)EMf->getBzTot(ix+1, iy+1, iz+1);
					  }

			 //Write VTK header
			 sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Magnetic Field from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d\n"
						   "VECTORS B float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (fieldtags[tagid].compare("E") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getEx(ix+1, iy+1, iz+1);
						  fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getEy(ix+1, iy+1, iz+1);
						  fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getEz(ix+1, iy+1, iz+1);
					  }

			 //Write VTK header
			 sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Electric Field from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS E float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (fieldtags[tagid].compare("Je") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 0);
						  fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 0);
						  fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 0);
					  }

			 //Write VTK header
			 sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Electron current from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS Je float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (fieldtags[tagid].compare("Ji") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  fieldwritebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 1);
						  fieldwritebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 1);
						  fieldwritebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 1);
					  }

			 //Write VTK header
			 sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Ion current from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS Ji float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);
		 }

		 if(EMf->isLittleEndian()){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  ByteSwap((unsigned char*) &fieldwritebuffer[iz][iy][ix][0],4);
						  ByteSwap((unsigned char*) &fieldwritebuffer[iz][iy][ix][1],4);
						  ByteSwap((unsigned char*) &fieldwritebuffer[iz][iy][ix][2],4);
					  }
		 }

		  int nelem = strlen(header);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_"<< fieldtags[tagid] << "_" << cycle << ".vtk";
		  MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &(fhArr[counter]));

		  if (vct->getCartesian_rank()==0){
			  MPI_Status   status;
			  MPI_File_write(fhArr[counter], header, nelem, MPI_BYTE, &status);
		  }

		  error_code = MPI_File_set_view(fhArr[counter], disp, EMf->getXYZeType(), EMf->getProcviewXYZ(), "native", MPI_INFO_NULL);
	      if(error_code!= MPI_SUCCESS){
              char error_string[100];
              int length_of_error_string, error_class;
              MPI_Error_class(error_code, &error_class);
              MPI_Error_string(error_class, error_string, &length_of_error_string);
              dprintf("Error in MPI_File_set_view: %s\n", error_string);
	      }

	      error_code = MPI_File_write_all_begin(fhArr[counter],fieldwritebuffer[0][0][0],(nxn-3)*(nyn-3)*(nzn-3),EMf->getXYZeType());
	      //error_code = MPI_File_iwrite(fhArr[counter], fieldwritebuffer[0][0][0],(nxn-3)*(nyn-3)*(nzn-3),EMf->getXYZeType(), &(requestArr[counter]));
	      if(error_code!= MPI_SUCCESS){
              char error_string[100];
              int length_of_error_string, error_class;
              MPI_Error_class(error_code, &error_class);
              MPI_Error_string(error_class, error_string, &length_of_error_string);
              dprintf("Error in MPI_File_iwrite: %s", error_string);
	      }
	      counter ++;
		}
	return counter;
}


int  WriteMomentsVTKNonblk(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct,int cycle,
		float*** momentswritebuffer,MPI_Request requestArr[14],MPI_File fhArr[14]){

	//All VTK output at grid cells excluding ghost cells
	const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
	const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
	const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
	const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
	const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
	const int    nPoints = dimX*dimY*dimZ;
	const string momentstags[]={"rho", "PXX", "PXY", "PXZ", "PYY", "PYZ", "PZZ"};
	const int    tagsize = size(momentstags);
	const string outputtag = col->getMomentsOutputTag();
	int counter=0,err;

	for(int tagid=0; tagid<tagsize; tagid++){
		 if (outputtag.find(momentstags[tagid], 0) == string::npos) continue;

		 for(int si=0;si<=1;si++){
			 momentswritebuffer = &(momentswritebuffer[counter]);
			 char  header[1024];
			 if (momentstags[tagid].compare("rho") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getRHOns(ix+1, iy+1, iz+1, si)*4*3.1415926535897;

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s density from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"rhoe":"rhoi");
			 }else if(momentstags[tagid].compare("PXX") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getpXXsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s pressure PXX from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PXXe":"PXXi");
			 }else if(momentstags[tagid].compare("PXY") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getpXYsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s pressure PXY from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PXYe":"PXYi");
			}else if(momentstags[tagid].compare("PXZ") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getpXZsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s pressure PXZ from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PXZe":"PXZi");
			}else if(momentstags[tagid].compare("PYY") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getpYYsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s pressure PYY from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PYYe":"PYYi");
			}else if(momentstags[tagid].compare("PYZ") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getpYZsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s pressure PYZ from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PYZe":"PYZi");
			}else if(momentstags[tagid].compare("PZZ") == 0){

				for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++)
						  momentswritebuffer[iz][iy][ix] = (float)EMf->getpZZsn(ix+1, iy+1, iz+1, si);

				//Write VTK header
				sprintf(header, "# vtk DataFile Version 2.0\n"
								   "%s pressure PZZ from iPIC3D\n"
								   "BINARY\n"
								   "DATASET STRUCTURED_POINTS\n"
								   "DIMENSIONS %d %d %d\n"
								   "ORIGIN 0 0 0\n"
								   "SPACING %f %f %f\n"
								   "POINT_DATA %d \n"
								   "SCALARS %s float\n"
								   "LOOKUP_TABLE default\n",(si%2==0)?"Electron":"Ion",dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,(si%2==0)?"PZZe":"PZZi");
			}

			 if(EMf->isLittleEndian()){
				 for(int iz=0;iz<nzn-3;iz++)
					  for(int iy=0;iy<nyn-3;iy++)
						  for(int ix= 0;ix<nxn-3;ix++){
							  ByteSwap((unsigned char*) &momentswritebuffer[iz][iy][ix],4);
						  }
			 }

		  int nelem = strlen(header);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_" << momentstags[tagid] << ((si%2==0)?"e":"i") << "_" << cycle << ".vtk";
		  MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhArr[counter]);

		  if (vct->getCartesian_rank()==0){
			  MPI_Status   status;
			  MPI_File_write(fhArr[counter], header, nelem, MPI_BYTE, &status);
		  }

	      int error_code = MPI_File_set_view(fhArr[counter], disp, MPI_FLOAT, EMf->getProcview(), "native", MPI_INFO_NULL);
	      if (error_code != MPI_SUCCESS) {
			char error_string[100];
			int length_of_error_string, error_class;
			MPI_Error_class(error_code, &error_class);
			MPI_Error_string(error_class, error_string, &length_of_error_string);
			dprintf("Error in MPI_File_set_view: %s\n", error_string);
		  }

	      error_code = MPI_File_write_all_begin(fhArr[counter],momentswritebuffer[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT);
	      //error_code = MPI_File_iwrite(fhArr[counter], momentswritebuffer[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &(requestArr[counter]));
	      if(error_code != MPI_SUCCESS){
	          char error_string[100];
	          int length_of_error_string, error_class;
	          MPI_Error_class(error_code, &error_class);
	          MPI_Error_string(error_class, error_string, &length_of_error_string);
	          dprintf("Error in MPI_File_iwrite: %s \n", error_string);
	      }
	      counter ++;
		 }//END OF SPECIES
	}
	return counter;
}



void ByteSwap(unsigned char * b, int n)
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void WriteTestPclsVTK(int nspec, Grid3DCU *grid, Particles3D *testpart, EMfields3D *EMf,
		CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,MPI_Request *testpartMPIReq, MPI_File *fh){
	/* the below is nonblocking collective IO
	 * const int nop = testpart[0].getNOP();

  	if(cycle>0){
  		MPI_Wait(headerReq, status);
  		MPI_Wait(dataReq, status);
  		MPI_Wait(footReq, status);
//  		MPI_File_close(&fh);
//  		int error_code=status->MPI_ERROR;
//  		if (error_code != MPI_SUCCESS) {
//  			char error_string[100];
//  			int length_of_error_string, error_class;
//
//  			MPI_Error_class(error_code, &error_class);
//  			MPI_Error_string(error_class, error_string, &length_of_error_string);
//  			dprintf("MPI_Wait error: %s\n", error_string);
//  		}
  	}else{
  		pclbuffersize = nop*3*1.2;
  		testpclPos = new float[pclbuffersize];
  	}

  	if(nop>pclbuffersize){
  		pclbuffersize = nop*3*1.2;
  		delete testpclPos;
  		testpclPos = new float[pclbuffersize];
  	}

//  	for (int pidx = 0; pidx < nop; pidx++) {
//  		testpclPos[pidx*3+0]=(float)(testpart[0]).getX(pidx);
//  		testpclPos[pidx*3+1]=(float)(testpart[0]).getY(pidx);
//  		testpclPos[pidx*3+2]=(float)(testpart[0]).getZ(pidx);
//  	}

  	//write to parallel vtk pvtu files
  	int  is=0;
	ostringstream filename;
	filename << col->getSaveDirName() << "/" << col->getSimName() << "_testparticle"<< testpart[is].get_species_num() << "_cycle" << cycle << ".vtu";
	MPI_File_open(vct->getComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

	  ofstream myfile;
	  myfile.open ("example.vtu");
	  myfile <<  "<?xml version=\"1.0\"?>\n"
				"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			    "  <UnstructuredGrid>\n"
				"    <Piece NumberOfPoints=\"1\" NumberOfCells=\"1\">\n"
				"		<Cells>\n"
				"			<DataArray type=\"UInt8\" Name=\"connectivity\" format=\"ascii\">0 1</DataArray>\n"
				"			<DataArray type=\"UInt8\" Name=\"offsets\" 		format=\"ascii\">1</DataArray>\n"
				"			<DataArray type=\"UInt8\" Name=\"types\"    	format=\"ascii\">1</DataArray>\n"
				"		</Cells>\n"
				"		<Points>\n"
				"        	<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"acscii\">\n" <<
				(testpart[0]).getX(100) << (testpart[0]).getY(100)  << (testpart[0]).getZ(100) <<
				 "			</DataArray>\n"
				  	  					"		</Points>\n"
				  	  					"	</Piece>\n"
				  	  					"	</UnstructuredGrid>\n"
				  	  					"</VTKFile>";
	  myfile.close();

  	char header[8192];
  	sprintf(header, "<?xml version=\"1.0\"?>\n"
  					"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n"
  				    "  <UnstructuredGrid>\n"
  					"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"1\">\n"
  					"		<Cells>\n"
  					"			<DataArray type=\"UInt8\" Name=\"connectivity\" format=\"ascii\">0 1</DataArray>\n"
  					"			<DataArray type=\"UInt8\" Name=\"offsets\" 		format=\"ascii\">1</DataArray>\n"
  					"			<DataArray type=\"UInt8\" Name=\"types\"    	format=\"ascii\">1</DataArray>\n"
  					"		</Cells>\n"
  					"		<Points>\n"
  					"        	<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n",
  					(EMf->isLittleEndian() ?"LittleEndian":"BigEndian"),3);

  	int nelem = strlen(header);
  	int charsize=sizeof(char);
  	MPI_Offset disp = nelem*charsize;

  	//MPI_File_iwrite(fh, header, nelem, MPI_BYTE, headerReq);
  	MPI_File_write(fh, header, nelem, MPI_BYTE, status);

  	int err = MPI_File_set_view(fh, disp, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
  	if(err){
  		          dprintf("Error in MPI_File_set_view\n");
  		      }

	//dprintf("testpart[is].getNOP() = %d, sizeof(SpeciesParticle)=%d, u = %f, x = %f ",testpart[0].getNOP(), sizeof(SpeciesParticle), pcl.get_u(), pcl.get_x());
	//const SpeciesParticle *pclptr = testpart[is].get_pclptr(0);
	//const SpeciesParticle * pclptr = (SpeciesParticle * )temppcl;
	//dprintf("temppcl u = %f, x= %f",pclptr->get_u() , pclptr->get_x());
//	MPI_Datatype particleType;
//	MPI_Type_vector (10,3,sizeof(SpeciesParticle),MPI_DOUBLE,&particleType);//testpart[0].getNOP()
//	MPI_Type_commit(&particleType);

  	//MPI_File_iwrite(fh, pclptr, 1, particleType, dataReq);
  	testpclPos[0]=1.1;testpclPos[1]=2.1;testpclPos[2]=3.1;
  	testpclPos[3]=1.1;testpclPos[4]=2.1;testpclPos[5]=3.1;
  	testpclPos[6]=1.1;testpclPos[7]=2.1;testpclPos[8]=3.1;
  	MPI_File_write_all(fh, testpclPos, 9, MPI_FLOAT, status);
  	 int tcount=0;
	  MPI_Get_count(status, MPI_FLOAT, &tcount);
	  dprintf(" wrote %i MPI_FLOAT",  tcount);
	int error_code=status->MPI_ERROR;
	if (error_code != MPI_SUCCESS) {
		char error_string[100];
		int length_of_error_string, error_class;

		MPI_Error_class(error_code, &error_class);
		MPI_Error_string(error_class, error_string, &length_of_error_string);
		dprintf("MPI_File_write error: %s\n", error_string);
	}

  	char foot[8192];
  	sprintf(foot, "			</DataArray>\n"
  	  					"		</Points>\n"
  	  					"	</Piece>\n"
  	  					"	</UnstructuredGrid>\n"
  	  					"</VTKFile>");
  	//nelem = strlen(foot);
  	//MPI_File_set_view(fh, disp+9, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
  	//MPI_File_iwrite(fh, foot, nelem, MPI_BYTE, footReq);
  	//MPI_File_write(fh, foot, nelem, MPI_BYTE, status);

  	if(cycle==LastCycle()){
  		MPI_Wait(headerReq, status);
  		MPI_Wait(dataReq, status);
  		MPI_Wait(footReq, status);
  	    MPI_File_close(&fh);
  	}
	 * */
	MPI_Status  *status;
	if(cycle>0){dprintf("Previous writing done");
		MPI_Wait(testpartMPIReq, status);dprintf("Previous writing done");
		int error_code=status->MPI_ERROR;
		if (error_code != MPI_SUCCESS) {
			char error_string[100];
			int length_of_error_string, error_class;

			MPI_Error_class(error_code, &error_class);
			MPI_Error_string(error_class, error_string, &length_of_error_string);
			dprintf("MPI_Wait error: %s\n", error_string);
		}
		else{
			if (vct->getCartesian_rank()==0) MPI_File_close(fh);
			dprintf("Previous writing done");
		}
	}

	//buffering if no full
	const int timesteps = 10;
	int nopStep[timesteps];



	//write to parallel vtk pvtu files

	int 		 is=0;

	char header[12048];
	sprintf(header, "<?xml version=\"1.0\"?>\n"
					"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n"
				    "  <UnstructuredGrid>\n"
					"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"1\">\n"
					"		<Cells>\n"
					"			<DataArray type=\"UInt8\" Name=\"connectivity\" format=\"ascii\">0 1</DataArray>\n"
					"			<DataArray type=\"UInt8\" Name=\"offsets\" 		format=\"ascii\">0</DataArray>\n"
					"			<DataArray type=\"UInt8\" Name=\"types\"    	format=\"ascii\">11</DataArray>\n"
					"		</Cells>\n"
					"		<Points>\n"
					"        	<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"
					"			</DataArray>\n"
					"		</Points>\n"
					"	</Piece>\n"
					"	</UnstructuredGrid>\n"
					"</VTKFile>", (EMf->isLittleEndian() ?"LittleEndian":"BigEndian"),testpart[is].getNOP());

	int nelem = strlen(header);
	int charsize=sizeof(char);
	MPI_Offset disp = nelem*charsize;

	ostringstream filename;
	filename << col->getSaveDirName() << "/" << col->getSimName() << "_testparticle"<< testpart[is].get_species_num() << "_cycle" << cycle << ".vtu";
	MPI_File_open(vct->getFieldComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, fh);

	MPI_File_set_view(*fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
	if (vct->getCartesian_rank()==0){

		MPI_File_iwrite(*fh, header, nelem, MPI_BYTE, testpartMPIReq);

	}

	/*
	 *
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="9" NumberOfCells="1">
        <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 1 2 3 4 5 6 7 8
        </DataArray>
		<DataArray type="Int32" Name="offsets" format="ascii">
		 0
		</DataArray>

		<DataArray type="UInt8" Name="types" format="ascii">
			1
		</DataArray>
        </Cells>

      <PointData Scalars="testpartID">
        <DataArray type="UInt16" Name="testpartID" format="ascii">
          0 1 2 3 4 5 6 7 8
        </DataArray>

        <DataArray type="Float32" Name="testpartVelocity" NumberOfComponents="3" format="ascii">
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
          -1 -1 -1
        </DataArray>
      </PointData>

      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0 0 0 1 0 0 2
          0 1 0 0 1 1 0 1 2
          0 2 0 0 2 1 0 2 2
        </DataArray>
      </Points>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
	 * */

}
