/***************************************************************************
 deel2vtk.cpp  -  Convert iPic3D hdf5 output to vtk format
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 Current version      : Feb 2015 by Jan Deca - convert particle cycle to vtk
 ************************************************************************** */


/* USAGE
 assuming proc and part hdf output from iPic3D
 ./hdf2vtk BEJRSTgg start end step
 assembles vtk file from cycle "start" to cycle "stop" using a certain "step"
 B = magnetic field
 E = electric field
 J = current e and i
 R = density e and i
 S = density e1,e2,i1,i2 (4 species)
 T = stress tensor
 gg = generates one file per part.hdf for each species. 
 
 example: ./hdf2vtk B 0 10 2 -> creates 6 vtk files for B.
*/




#include <mpi.h>
#include "hdf5.h"
#include "Alloc.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}


// Functions
int readsettings();
void readvect(int cycle, string dir, string vectname, string speciesname, double ***VX, double ***EY,double ***EZ);
void readscalar(int cycle, string dir, string scalarname, string speciesname, double ***RHO);
void writevect(int cycle, string vectname, string speciesname, double ***EX, double ***EY,double ***EZ);
void writescalar(int cycle, string scalarname, string speciesname, double*** PHI);
void writedensity(int cycle, string scalarname, double ***RHOe, double ***RHOi);
void writedensity4(int cycle, string scalarname, double ***RHOe, double ***RHOi, double ***RHOe2, double ***RHOi2);
void writetensor(int cycle, string tensorname, string speciesname, double*** EXX, double*** EXY, double*** EXZ, double*** EYY, double*** EYZ, double*** EZZ);
void extract_pressure(double qom, double*** VX, double*** VY, double*** VZ, double*** N, double*** pXX, double*** pXY, double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ);
void readwriteparticles(int cycle, string dir, string speciesname);
void readwritepartproc(int cycle, string dir, string speciesname);

// Variables
int nxn, nyn, nzn;
int nxc, nyc, nzc;
double dx, dy, dz;
int XLEN, YLEN, ZLEN;

double Lx, Ly, Lz;
int ns;
double* qom;

int nproc;
double *temp_storageX;
double *temp_storageY;
double *temp_storageZ;
double *temp_storage;

hid_t	file_id;
hid_t	dataset_id;
herr_t	status;

double *temp_q;
double *temp_u;
double *temp_v;
double *temp_w;
double *temp_x;
double *temp_y;
double *temp_z;





// Main
int main (int argc, char **argv) {

	
	
	int start, end, step;
	char field [20];
	int startframe, endframe, nsteps;
	
	sscanf(argv[1],"%s",field);
	sscanf(argv[2],"%d",&start);
	sscanf(argv[3],"%d",&end);
	sscanf(argv[4],"%d",&step);

	string fields = string(field);
	cout << "Converting: " << fields << endl;
	//set up MPI
	int myrank, nprocs;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	nsteps = int(((end-start)/step + 1 + myrank)/nprocs);
	startframe = start;
	for (int i = 0; i<myrank; i++) {
		startframe = startframe + int(((end-start)/step + 1 + i)/nprocs)*step;
	}
	endframe = startframe + (nsteps-1)*step;
	
	
	// read "settings.hdf"
	int out;
	out = readsettings();
	if(out<0){
		return -1;
	}
	
	for (int ii=startframe;ii <=endframe; ii+=step){
		
		nxn = nxc/XLEN;
		nyn = nyc/YLEN;
		nzn = nzc/ZLEN;
		
		dx = Lx/nxc;
		dy = Ly/nyc;
		dz = Lz/nzc;
		
		
		temp_storageX = new double[(nxn+1)*(nyn+1)*(nzn+1)];
		temp_storageY = new double[(nxn+1)*(nyn+1)*(nzn+1)];
		temp_storageZ = new double[(nxn+1)*(nyn+1)*(nzn+1)];
		temp_storage = new double[(nxn+1)*(nyn+1)*(nzn+1)];
	
		double*** VX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** VY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** VZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TXX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TXY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TXZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TYY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TYZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TZZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** NE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** NI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
        double*** NE2 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
        double*** NI2 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
        double*** RHO = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);

		
		//E
		if(fields.find("E")!=string::npos){
			readvect(ii, "/fields/","E","", VX, VY, VZ);
			writevect(ii, "E", "", VX, VY, VZ);
		}
		//B
		if(fields.find("B")!=string::npos){
			readvect(ii,"/fields/","B","", VX, VY, VZ);
			writevect(ii, "B", "", VX, VY, VZ);
		}
		//Rho
		if(fields.find("R")!=string::npos || fields.find("T")!=string::npos){
			readscalar(ii,"/moments/species_0/","rho","e",  NE);
			readscalar(ii,"/moments/species_1/","rho","i",  NI);
			if(fields.find("R")!=string::npos){
				writedensity(ii,"rho", NE, NI);
			}
		}
        if(fields.find("S")!=string::npos){
            readscalar(ii,"/moments/species_0/","rho","e1",  NE);
            readscalar(ii,"/moments/species_1/","rho","i1",  NI);
            readscalar(ii,"/moments/species_2/","rho","e2",  NE2);
            readscalar(ii,"/moments/species_3/","rho","i2",  NI2);
            if(fields.find("S")!=string::npos){
                writedensity4(ii,"rho", NE, NI,NE2,NI2);
            }
        }

		//Je
		if(fields.find("J")!=string::npos || fields.find("T")!=string::npos){
			readvect(ii,"/moments/species_0/","J","e",  VX, VY, VZ);
			if(fields.find("J")!=string::npos){
				writevect(ii, "J", "e", VX, VY, VZ);
			}
		}
		//Pe
		if(fields.find("T")!=string::npos){
			readscalar(ii,"/moments/species_0/","pXX","e",  TXX);
			readscalar(ii,"/moments/species_0/","pXY","e",  TXY);
			readscalar(ii,"/moments/species_0/","pXZ","e",  TXZ);
			readscalar(ii,"/moments/species_0/","pYY","e",  TYY);
			readscalar(ii,"/moments/species_0/","pYZ","e",  TYZ);
			readscalar(ii,"/moments/species_0/","pZZ","e",  TZZ);
			extract_pressure(qom[0], VX, VY, VZ, NE, TXX, TXY, TXZ, TYY, TYZ, TZZ);
			writetensor(ii, "P", "e", TXX, TXY, TXZ, TYY, TYZ, TZZ);
		}
		//Ji
		if(fields.find("J")!=string::npos || fields.find("T")!=string::npos){
			readvect(ii,"/moments/species_1/","J","i",  VX, VY, VZ);
			if(fields.find("J")!=string::npos){
				writevect(ii, "J", "i", VX, VY, VZ);
			}
		}
		//Pi
		if(fields.find("T")!=string::npos){
			readscalar(ii,"/moments/species_1/","pXX","i",  TXX);
			readscalar(ii,"/moments/species_1/","pXY","i",  TXY);
			readscalar(ii,"/moments/species_1/","pXZ","i",  TXZ);
			readscalar(ii,"/moments/species_1/","pYY","i",  TYY);
			readscalar(ii,"/moments/species_1/","pYZ","i",  TYZ);
			readscalar(ii,"/moments/species_1/","pZZ","i",  TZZ);
			extract_pressure(qom[1], VX, VY, VZ, NI, TXX, TXY, TXZ, TYY, TYZ, TZZ);
			writetensor(ii, "P", "i", TXX, TXY, TXZ, TYY, TYZ, TZZ);		
		}
		//PHI
		if(fields.find("F")!=string::npos){
			readscalar(ii,"/potentials/","phi","",  NE);
		}

		//Particles
		if(fields.find("dd0")!=string::npos){
		  cout << "we gonna read particles" << endl;
		  readwriteparticles(ii,"/particles/species_0","electrons");
		}
		if(fields.find("dd1")!=string::npos){
		  cout << "we gonna read particles" << endl;
		  readwriteparticles(ii,"/particles/species_1","ions");
		}
		if(fields.find("gg")!=string::npos){
		  readwritepartproc(ii,"/particles/species_0","e");
		  readwritepartproc(ii,"/particles/species_1","i");

		}
		
		
		if(ii==endframe){
			delArr3(VX,nxn*XLEN,nyn*YLEN);
			delArr3(VY,nxn*XLEN,nyn*YLEN);
			delArr3(VZ,nxn*XLEN,nyn*YLEN);
			delArr3(TXX,nxn*XLEN,nyn*YLEN);
			delArr3(TXY,nxn*XLEN,nyn*YLEN);
			delArr3(TXZ,nxn*XLEN,nyn*YLEN);
			delArr3(TYY,nxn*XLEN,nyn*YLEN);
			delArr3(TYZ,nxn*XLEN,nyn*YLEN);
			delArr3(TZZ,nxn*XLEN,nyn*YLEN);
            delArr3(NE,nxn*XLEN,nyn*YLEN);
            delArr3(NI,nxn*XLEN,nyn*YLEN);
            delArr3(NE2,nxn*XLEN,nyn*YLEN);
            delArr3(NI2,nxn*XLEN,nyn*YLEN);
            delArr3(RHO,nxn*XLEN,nyn*YLEN);
			delete[] temp_storageX;
			delete[] temp_storageY;
			delete[] temp_storageZ;
			delete[] temp_storage;
            delete[] temp_q;
            delete[] temp_u;
            delete[] temp_v;
            delete[] temp_w;
            delete[] temp_x;
            delete[] temp_y;
            delete[] temp_z;

		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return(0);
}



// read "settings.hdf"
int readsettings(){
	
	// Open the  settings file
	file_id = H5Fopen("settings.hdf", H5F_ACC_RDWR, H5P_DEFAULT);
	
	if (file_id < 0){
	
		cout << "couldn't open file: settings.hdf" << endl;
		return -1;
	} 
	else {
		// First read the topology
		int nproc;
		dataset_id = H5Dopen1(file_id, "/topology/Nprocs");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nproc);
		status = H5Dclose(dataset_id);
		
		dataset_id = H5Dopen1(file_id, "/topology/XLEN");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&XLEN);
		status = H5Dclose(dataset_id);
		
		dataset_id = H5Dopen1(file_id, "/topology/YLEN");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&YLEN);
		status = H5Dclose(dataset_id);
		
		dataset_id = H5Dopen1(file_id, "/topology/ZLEN");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ZLEN);
		status = H5Dclose(dataset_id);
		
		// read Lx
		dataset_id = H5Dopen1(file_id, "/collective/Lx");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
		status = H5Dclose(dataset_id);
		// read Ly
		dataset_id = H5Dopen1(file_id, "/collective/Ly");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
		status = H5Dclose(dataset_id);
		// read Lz
		dataset_id = H5Dopen1(file_id, "/collective/Lz");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lz);
		status = H5Dclose(dataset_id);
		// read nxc
		dataset_id = H5Dopen1(file_id, "/collective/Nxc");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nxc);
		status = H5Dclose(dataset_id);
		// read nyc
		dataset_id = H5Dopen1(file_id, "/collective/Nyc");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nyc);
		status = H5Dclose(dataset_id);
		// read nyc
		dataset_id = H5Dopen1(file_id, "/collective/Nzc");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nzc);
		status = H5Dclose(dataset_id);
		// read ns
		dataset_id = H5Dopen1(file_id, "/collective/Ns");
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ns);
		// read qom
		qom = new double[ns];
		stringstream specie;
		string temp;
		for (int is=0; is<ns; is++){
			specie.clear();
			specie.str("");
			specie << is;
			temp = "/collective/species_"+specie.str()+"/qom";
			dataset_id = H5Dopen1(file_id, temp.c_str());
			status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&qom[is]);
		}
		
		// close settings
		status = H5Fclose(file_id);
		return 0;
	}
}
// read vector
void readvect(int cycle, string dir, string vectname, string speciesname, double*** EX, double*** EY,double*** EZ) {
	
	hid_t proc_file_id;
	stringstream cc;
	cout << "READING VECTOR " << vectname + speciesname<< " FROM HDF5 FILES - Time Level = "<< cycle << endl;
	
    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() << endl;
				temp = "proc" + ss.str() + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = dir+vectname+"x/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen1(proc_file_id,temp.c_str());
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
				temp = dir+vectname+"y/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen1(proc_file_id,temp.c_str());
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
				temp = dir+vectname+"z/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen1(proc_file_id,temp.c_str());
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
				int node=0;
				for (int ii=0; ii < (nxn+1);ii++)
					for (int jj=0; jj < (nyn+1);jj++)
						for (int kk=0; kk < (nzn+1);kk++){
							if (ii!=nxn && jj!= nyn && kk!= nzn){
								EX[ii + nxn*i][jj + nyn*j][kk + nzn*k] = temp_storageX[node];
								EY[ii + nxn*i][jj + nyn*j][kk + nzn*k] = temp_storageY[node];
								EZ[ii + nxn*i][jj + nyn*j][kk + nzn*k] = temp_storageZ[node];
							}
							node++;
						}
				H5Fclose(proc_file_id);
			}
}
// read scalar
void readscalar(int cycle, string dir, string scalarname, string speciesname, double*** RHO) {
	
	hid_t proc_file_id;
	stringstream cc;
	cout << "READING SCALAR " << scalarname+speciesname << "  FROM HDF5 FILES - Time Level= "<< cycle << endl;
	
    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
//	int i=XLEN/2;
//	int k=ZLEN/2;
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() << endl;
				temp = "proc" + ss.str() + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = dir+scalarname+"/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen1(proc_file_id,temp.c_str());
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storage); status = H5Dclose(dataset_id);
				
				int node=0;
				for (int ii=0; ii < (nxn+1);ii++)
					for (int jj=0; jj < (nyn+1);jj++)
						for (int kk=0; kk < (nzn+1);kk++){
							if (ii!=nxn && jj!= nyn && kk!= nzn){
								RHO[ii + nxn*i][jj + nyn*j][kk + nzn*k] = temp_storage[node];
							}
							node++;
						}
				H5Fclose(proc_file_id);
			}
}
// write vector
void writevect(int cycle, string vectname, string speciesname, double*** EX, double*** EY,double*** EZ) {
	
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +speciesname +"_Jan_cycle" +stringcycle.str();
	temp += ".vtk";
	//cout << "Writing file: " << temp << endl;
	ofstream my_file(temp.c_str());
	//cout << "writing to file mesh points for" << vectname+addname << endl;
	my_file << "# vtk DataFile Version 1.0" << endl;
	my_file << vectname << " Field from Parsek" << endl;
	my_file << "ASCII" << endl;
	my_file << "DATASET STRUCTURED_POINTS" << endl;
	my_file << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_file << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_file << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_file << "VECTORS " << vectname+speciesname << " float" << endl;
	cout << "WRITING VECTOR " << vectname+speciesname<<" TO VTK FILE - Time Level= "<< cycle  << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EX[ii][jj][kk] << " " << EY[ii][jj][kk] << " " << EZ[ii][jj][kk] << endl;
			}
	my_file.close();
}
// write scalar
void writescalar(int cycle, string scalarname, string speciesname, double*** PHI) {
	
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname + speciesname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
//	cout << "Writing file: " << temp << endl;
	ofstream my_file(temp.c_str());
//	cout << "writing to file mesh points for " << scalarname+addname << endl;
	my_file << "# vtk DataFile Version 1.0" << endl;
	my_file << scalarname << " Field from Parsek" << endl;
	my_file << "ASCII" << endl;
	my_file << "DATASET STRUCTURED_POINTS" << endl;
	my_file << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_file << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_file << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_file << "SCALARS " << scalarname+speciesname << " float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname +speciesname<<" TO VTK FILE - Time Level= "<< cycle  << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << PHI[ii][jj][kk] << endl;
			}
	my_file.close();
}
//write density
void writedensity(int cycle, string scalarname, double*** RHOe, double*** RHOi) {
	
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_Jan_cycle" +stringcycle.str();
	temp += ".vtk";
//	cout << "Writing file: " << temp << endl;
	ofstream my_file(temp.c_str());
//	cout << "writing to file mesh points for " << scalarname << endl;
	my_file << "# vtk DataFile Version 1.0" << endl;
	my_file << scalarname << " Field from Parsek" << endl;
	my_file << "ASCII" << endl;
	my_file << "DATASET STRUCTURED_POINTS" << endl;
	my_file << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_file << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_file << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_file << "SCALARS " << scalarname << "0 float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"e TO VTK FILE - Time Level= "<< cycle  << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << RHOe[ii][jj][kk] << endl;
			}
	my_file << "SCALARS " << scalarname << "1 float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"i TO VTK FILE - Time Level= "<< cycle  << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << RHOi[ii][jj][kk] << endl;
			}
	my_file.close();
}
//write tensor
void writetensor(int cycle, string tensorname, string speciesname, double*** EXX, double*** EXY, double*** EXZ, double*** EYY, double*** EYZ, double*** EZZ) {

	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = tensorname + speciesname +"_Jan_cycle" +stringcycle.str();
	temp += ".vtk";
//	cout << "Writing file: " << temp << endl;
	ofstream my_file(temp.c_str());
//	cout << "writing to file mesh points for " << tensorname +addname << endl;
	my_file << "# vtk DataFile Version 1.0" << endl;
	my_file << tensorname +speciesname << " Tensor from Parsek" << endl;
	my_file << "ASCII" << endl;
	my_file << "DATASET STRUCTURED_POINTS" << endl;
	my_file << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_file << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_file << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	
    cout << "WRITING Tensor " << tensorname + speciesname<<" TO VTK FILE - Time Level= "<< cycle  << endl;
	
	my_file << "SCALARS " << tensorname+ speciesname << "xx float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EXX[ii][jj][kk] << endl;
			}
	my_file << "SCALARS " << tensorname+speciesname << "xy float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EXY[ii][jj][kk] << endl;
			}
	my_file << "SCALARS " << tensorname+speciesname << "xz float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EXZ[ii][jj][kk] << endl;
			}
	my_file << "SCALARS " << tensorname+speciesname << "yy float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EYY[ii][jj][kk] << endl;
			}
	my_file << "SCALARS " << tensorname+speciesname << "yz float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EYZ[ii][jj][kk] << endl;
			}
	my_file << "SCALARS " << tensorname+speciesname << "zz float" << endl;
	my_file << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_file << EZZ[ii][jj][kk] << endl;
			}
	my_file.close();
}
//correct tensor for flow velocity
void extract_pressure(double qom, double*** VX, double*** VY, double*** VZ, double*** N, double*** pXX, double*** pXY, double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ) {
	
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				if(N[ii][jj][kk]!=0.0)
					pXX[ii][jj][kk] = pXX[ii][jj][kk]-VX[ii][jj][kk]*VX[ii][jj][kk]/N[ii][jj][kk];
				pXX[ii][jj][kk] = pXX[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pXY[ii][jj][kk] = pXY[ii][jj][kk]-VX[ii][jj][kk]*VY[ii][jj][kk]/N[ii][jj][kk];
				pXY[ii][jj][kk] = pXY[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pXZ[ii][jj][kk] = pXZ[ii][jj][kk]-VX[ii][jj][kk]*VZ[ii][jj][kk]/N[ii][jj][kk];
				pXZ[ii][jj][kk] = pXZ[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pYY[ii][jj][kk] = pYY[ii][jj][kk]-VY[ii][jj][kk]*VY[ii][jj][kk]/N[ii][jj][kk];
				pYY[ii][jj][kk] = pYY[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pYZ[ii][jj][kk] = pYZ[ii][jj][kk]-VY[ii][jj][kk]*VZ[ii][jj][kk]/N[ii][jj][kk];
				pYZ[ii][jj][kk] = pYZ[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pZZ[ii][jj][kk] = pZZ[ii][jj][kk]-VZ[ii][jj][kk]*VZ[ii][jj][kk]/N[ii][jj][kk];
				pZZ[ii][jj][kk] = pZZ[ii][jj][kk] / qom;
			}
}


void readwriteparticles(int cycle, string dir, string speciesname) {

  hid_t part_file_id, dataset_id;
  hid_t dspace_id;
  stringstream cc;
  cout << "READING " << speciesname << " FROM HDF5 FILES - Time level = " << cycle << endl;

  int numpartproc[XLEN*YLEN*ZLEN];

  
  cc.clear();
  cc.str("");
  cc << cycle;
  string temp;
  int part;
  for (int i=0; i<XLEN; i++)
    for (int j=0; j<YLEN; j++)
      for (int k=0; k<ZLEN; k++) {
	//cout << "i=" << i << " j=" << j << " k=" << k << endl;
	part = i*YLEN*ZLEN + j*ZLEN + k;
	stringstream ss;
	ss << part;
	//cout << "ss=" << ss.str() << endl;
	temp = "part" + ss.str() + ".hdf";
	part_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	//read data
		temp = dir + "/q/cycle_" + cc.str();
	//	cout << "dataset = " << temp << endl;
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	//	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
	//cout << H5Dget_space(dataset_id) << endl;
	dspace_id = H5Dget_space(dataset_id);
	const int ndims = H5Sget_simple_extent_ndims(dspace_id);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dspace_id, dims, NULL);
	cout << dims[0] << endl;
	numpartproc[part] = dims[0];

	H5Fclose(part_file_id);
      }
  //cout << numpartproc[0] << numpartproc[1] << numpartproc[2] << numpartproc[3] << endl;
  int totpart = 0;
  for (int a=0; a<XLEN*YLEN*ZLEN; a++) {
    totpart+=numpartproc[a];
  }

  stringstream stringcycle;
  stringcycle << cycle;
  string opfn;
  opfn = "particles_" + speciesname + "_Jan_cycle" + stringcycle.str() + ".vtk";
  ofstream my_file(opfn.c_str());
  
  cc.clear();
  cc.str("");
  cc << cycle;
  for (int i=0; i<XLEN; i++)
    for (int j=0; j<YLEN; j++)
      for (int k=0; k<ZLEN; k++) {
	part = i*YLEN*ZLEN + j*ZLEN + k;

	temp_q = new double[numpartproc[part]];
	temp_u = new double[numpartproc[part]];
	temp_v = new double[numpartproc[part]];
	temp_w = new double[numpartproc[part]];
	temp_x = new double[numpartproc[part]];
	temp_y = new double[numpartproc[part]];
	temp_z = new double[numpartproc[part]];
	cout << "READING part" << part << ".hdf" << endl;
	stringstream ss;
	ss << part;
	temp = "part" + ss.str() + ".hdf";
	part_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	temp = dir + "/q/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_q);
	status = H5Dclose(dataset_id);
	temp = dir + "/u/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_u);
	status = H5Dclose(dataset_id);
	temp = dir + "/v/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_v);
	status = H5Dclose(dataset_id);
	temp = dir + "/w/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_w);
	status = H5Dclose(dataset_id);
	temp = dir + "/x/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_x);
	status = H5Dclose(dataset_id);
	temp = dir + "/y/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_y);
	status = H5Dclose(dataset_id);
	temp = dir + "/z/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_z);
	status = H5Dclose(dataset_id);
	
	cout << "WRITING part" << part << ".hdf" << endl;
	for (int aa=0; aa< numpartproc[part]; aa++) {
	  my_file << temp_q[aa] << " " << temp_x[aa] << " " << temp_y[aa] << " " << temp_z [aa] << " " << temp_u[aa] << " " << temp_v[aa] << " " << temp_w[aa] << endl;
	}


	
	H5Fclose(part_file_id);
      }

  my_file.close();
  
}

void readwritepartproc(int cycle, string dir, string speciesname) {

  hid_t part_file_id, dataset_id;
  hid_t dspace_id;
  stringstream cc;
  if (patch::to_string(speciesname) == "e") {
    cout << "READING electrons FROM HDF5 FILES - Time level = " << cycle << endl;
  } else if (patch::to_string(speciesname) == "i") {
    cout << "READING ions FROM HDF5 FILES - Time level = " << cycle << endl;
  } else {
    cout << "Ai Caramba, som-a-think-a is-a wronk!" << endl;
  }

  int numpartproc; //[XLEN*YLEN*ZLEN];

  
  cc.clear();
  cc.str("");
  cc << cycle;
  string temp;
  int part;
  for (int i=0; i<XLEN; i++)
    for (int j=0; j<YLEN; j++)
      for (int k=0; k<ZLEN; k++) {
	//cout << "i=" << i << " j=" << j << " k=" << k << endl;
	part = i*YLEN*ZLEN + j*ZLEN + k;
	stringstream ss;
	ss << part;
	//cout << "ss=" << ss.str() << endl;
	temp = "part" + ss.str() + ".hdf";
	part_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	//read data
		temp = dir + "/q/cycle_" + cc.str();
	//	cout << "dataset = " << temp << endl;
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	//	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
	//cout << H5Dget_space(dataset_id) << endl;
	dspace_id = H5Dget_space(dataset_id);
	const int ndims = H5Sget_simple_extent_ndims(dspace_id);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dspace_id, dims, NULL);
	//	cout << dims[0] << endl;
	numpartproc = dims[0];
    status = H5Dclose(dataset_id);  //just added
          
	temp_q = new double[numpartproc];
	temp_u = new double[numpartproc];
	temp_v = new double[numpartproc];
	temp_w = new double[numpartproc];
	temp_x = new double[numpartproc];
	temp_y = new double[numpartproc];
	temp_z = new double[numpartproc];
	cout << "READING part" << part << ".hdf" << endl;
	//	stringstream ss;
	//ss << part;
	temp = "part" + ss.str() + ".hdf";
	part_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	temp = dir + "/q/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_q);
	status = H5Dclose(dataset_id);
	temp = dir + "/u/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_u);
	status = H5Dclose(dataset_id);
	temp = dir + "/v/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_v);
	status = H5Dclose(dataset_id);
	temp = dir + "/w/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_w);
	status = H5Dclose(dataset_id);
	temp = dir + "/x/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_x);
	status = H5Dclose(dataset_id);
	temp = dir + "/y/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_y);
	status = H5Dclose(dataset_id);
	temp = dir + "/z/cycle_" + cc.str();
	dataset_id = H5Dopen1(part_file_id,temp.c_str());
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_z);

	H5Fclose(part_file_id);

	stringstream stringcycle;
	stringcycle << cycle;
	string opfn;
	
	
	opfn = "part_" + speciesname + "_Jan_cycle" + stringcycle.str() + "_" + patch::to_string(i) + "_" + patch::to_string(j) + "_" + patch::to_string(k) + "_" + patch::to_string(part) + ".vtp";
	ofstream my_file(opfn.c_str());

	cout << "WRITING " << opfn << endl;
	my_file << speciesname << " q x y z u v w " << endl; 
	my_file << "proctopology " << patch::to_string(part) << " " << patch::to_string(XLEN) << " " << patch::to_string(YLEN) << " " << patch::to_string(ZLEN) << " " << patch::to_string(i) << " " << patch::to_string(j) << " " << patch::to_string(k) << endl;
	for (int aa=0; aa< numpartproc; aa++) {
	  my_file << temp_q[aa] << " " << temp_x[aa] << " " << temp_y[aa] << " " << temp_z [aa] << " " << temp_u[aa] << " " << temp_v[aa] << " " << temp_w[aa] << endl;
	}

	my_file.close();
      }
}

//write density - 4 species
void writedensity4(int cycle, string scalarname, double*** RHOe, double*** RHOi, double*** RHOe2, double*** RHOi2) {
    string temp;
    stringstream stringcycle;
    stringcycle << cycle;
    temp = scalarname +"4S_Jan_cycle" +stringcycle.str();
    temp += ".vtk";
    //	cout << "Writing file: " << temp << endl;
    ofstream my_file(temp.c_str());
    //	cout << "writing to file mesh points for " << scalarname << endl;
    my_file << "# vtk DataFile Version 1.0" << endl;
    my_file << scalarname << " Field from Parsek" << endl;
    my_file << "ASCII" << endl;
    my_file << "DATASET STRUCTURED_POINTS" << endl;
    my_file << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
    my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    my_file << "SPACING " << dy << " " << dy << " " << dz << endl;
    my_file << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
    my_file << "SCALARS " << scalarname << "e1 float" << endl;
    my_file << "LOOKUP_TABLE default" << endl;
    cout << "WRITING SCALAR " << scalarname <<"e1 TO VTK FILE - Time Level= "<< cycle  << endl;
    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){
                my_file << RHOe[ii][jj][kk] << endl;
            }
    my_file << "SCALARS " << scalarname << "i1 float" << endl;
    my_file << "LOOKUP_TABLE default" << endl;
    cout << "WRITING SCALAR " << scalarname <<"i1 TO VTK FILE - Time Level= "<< cycle  << endl;
    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){
                my_file << RHOi[ii][jj][kk] << endl;
            }
    my_file << "SCALARS " << scalarname << "e2 float" << endl;
    my_file << "LOOKUP_TABLE default" << endl;
    cout << "WRITING SCALAR " << scalarname <<"e2 TO VTK FILE - Time Level= "<< cycle  << endl;
    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){
                my_file << RHOe2[ii][jj][kk] << endl;
            }
    my_file << "SCALARS " << scalarname << "i2 float" << endl;
    my_file << "LOOKUP_TABLE default" << endl;
    cout << "WRITING SCALAR " << scalarname <<"i2 TO VTK FILE - Time Level= "<< cycle  << endl;
    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){
                my_file << RHOi2[ii][jj][kk] << endl;
            }
    my_file.close();
}





















