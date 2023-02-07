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


#include "mpi.h"
#include "MPIdata.h"
#include "iPic3D.h"
#include "TimeTasks.h"
#include "ipicdefs.h"
#include "debug.h"
#include "Parameters.h"
#include "ompdefs.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"
#include "Particles3Dcomm.h"
#include "Timing.h"
#include "ParallelIO.h"
#include "Collisions.h"
//
#ifndef NO_HDF5
#include "OutputWrapperFPP.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>

#include "Moments.h" // for debugging

#ifdef USE_CATALYST_LEGACY
#include "Adaptor_legacy.h"
#endif

using namespace iPic3D;
//MPIdata* iPic3D::c_Solver::mpi=0;

c_Solver::~c_Solver()
{
  delete col; // configuration parameters ("collectiveIO")
  delete vct; // process topology
  delete grid; // grid
  delete EMf; // field
  delete colls; // Collisions
#ifndef NO_HDF5
  delete outputWrapperFPP;
#endif
  // delete particles
  //
  if(part)
  {
    for (int i = 0; i < ns; i++)
    {
      // placement delete
      part[i].~Particles3D();
    }
    free(part);
  }

  #ifdef USE_CATALYST_LEGACY
  Adaptor_legacy::Finalize();
  #endif

  delete [] Ke;
  delete [] rho;
  delete [] momentum;
  delete [] Count;
  delete [] Qdel;
  delete [] Qrep;
  delete [] Qexo;
  delete my_clock;
}


/*  -------------- */
/*!  Initialization  */
/*  -------------- */
int c_Solver::Init(int argc, char **argv) {
  #if defined(__MIC__)
  assert_eq(DVECWIDTH,8);
  #endif

  Parameters::init_parameters();

  nprocs = MPIdata::get_nprocs();
  myrank = MPIdata::get_rank();

  col = new Collective(argc, argv);

  verbosity = (col->getVerbose())*(myrank==0);

  if (verbosity)
    dprintf("(1) New object Collective created");
  
  ns = col->getNs();            
  restart_cycle  = col->getRestartOutputCycle();
  SaveDirName    = col->getSaveDirName();
  RestartDirName = col->getRestartDirName();
  restart_status = col->getRestart_status();
  first_cycle = col->getLast_cycle()+1;

  vct = new VCtopology3D(*col);

  if (verbosity)
    dprintf("(2) New object Virtual Cartesian topology created");

  if (nprocs != vct->getNprocs() and myrank==0)
    eprintf("Error: %i processes cant be mapped into %ix%ix%i matrix: Change XLEN,YLEN, ZLEN in method VCtopology3D.init()", nprocs, vct->getXLEN(), vct->getYLEN(), vct->getZLEN())

  vct->setup_vctopology(MPIdata::get_PicGlobalComm());

  if (myrank == 0)
  {
    MPIdata::instance().Print();
    vct->Print();
    col->save();
    if (verbosity)	   
      col->Print();
  }

  grid = new Grid3DCU(col, vct);

  if (verbosity)
    dprintf("(3) Done with local grid. Make EM field object.");

  EMf = new EMfields3D(col, grid, vct);
 
  if (col->getcollisionProcesses())
  {
    colls = new Collisions(col, vct, grid, EMf);
    if (verbosity)
      dprintf("(3a) create new collison object.");
  }

  EMf->initDipole();

  if (verbosity)
    dprintf("(4) EMF initialized with initDipole() method. Beginning allocation of particles.");

  part = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int i = 0; i < ns; i++)
  {
    new(&part[i]) Particles3D(i,col,vct,grid);
    if (restart_status == 0) 
    {
      part[i].maxwellianDipole(EMf,col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
      part[i].reserve_remaining_particle_IDs();
    }
    else
    {
      part[i].load_restart_pcls();
    }
  }

  if (verbosity)
    dprintf("(5) Particle initialization Dipole or loading from restart completed.");

  //allocate test particles if any
  nstestpart = col->getNsTestPart();
  if(nstestpart>0 && restart_status==0)
  {
    testpart = (Particles3D*) malloc(sizeof(Particles3D)*nstestpart);
    for (int i = 0; i < nstestpart; i++)
    {
      new(&testpart[i]) Particles3D(i+ns,col,vct,grid);//species id for test particles is increased by ns
      testpart[i].pitch_angle_energy(EMf);
    }
  }
  // TBR this part above on test pcls --> careful to init_output_files

  #ifndef NO_HDF5
  if (col->getWriteMethod() == "shdf5" || restart_cycle>0 || (col->getWriteMethod()=="pvtk" && !col->particle_output_is_off()) )
    {
      outputWrapperFPP = new OutputWrapperFPP;
      fetch_outputWrapperFPP().init_output_files(col,vct,grid,EMf,part,ns,testpart,nstestpart);
 
      if (verbosity)
        dprintf("(5a) Output file. Created restart*.hdf (empty) and setting.hdf (full).");   
    }
  #endif
 
  if(!col->field_output_is_off())
  {
    if(col->getWriteMethod()=="pvtk")
    {
      if(!(col->getFieldOutputTag()).empty())
        fieldwritebuffer = newArr4(float,(grid->getNZN()-3),grid->getNYN()-3,grid->getNXN()-3,3);
      if(!(col->getMomentsOutputTag()).empty())
        momentwritebuffer=newArr3(float,(grid->getNZN()-3), grid->getNYN()-3, grid->getNXN()-3);
    }
    else if(col->getWriteMethod()=="nbcvtk")
    {
      momentreqcounter=0;
      fieldreqcounter = 0;
      if(!(col->getFieldOutputTag()).empty())
        fieldwritebuffer = newArr4(float,(grid->getNZN()-3)*4,grid->getNYN()-3,grid->getNXN()-3,3);
      if(!(col->getMomentsOutputTag()).empty())
        momentwritebuffer=newArr3(float,(grid->getNZN()-3)*14, grid->getNYN()-3, grid->getNXN()-3);
    }
  }
	  
  if(!col->spectra_output_is_off())
  {
    if(col->getWriteMethod()=="pvtk")
    {
      if(!(col->getSpectraOutputTag()).empty())
      {
        if ((grid->getNXN()-3)%col->getDeltaX() != 0 || (grid->getNYN()-3)%col->getDeltaY() != 0 || (grid->getNZN()-3)%col->getDeltaZ() != 0)
          eprintf("ERROR: Number of cells in MPI subdomain not a multiple of Delta (Energy Spectra).")
        int Nspece = (col->getEende() - col->getEstarte())/col->getdEe() + 1.001;
        int Nspeci = (col->getEendi() - col->getEstarti())/col->getdEi() + 1.001;
        spectrawritebuffere=newArr4(float,(grid->getNZN()-3)/(col->getDeltaZ()),(grid->getNYN()-3)/(col->getDeltaY()),(grid->getNXN()-3)/(col->getDeltaX()),Nspece+1);
        spectrawritebufferi=newArr4(float,(grid->getNZN()-3)/(col->getDeltaZ()),(grid->getNYN()-3)/(col->getDeltaY()),(grid->getNXN()-3)/(col->getDeltaX()),Nspeci+1);
      }
    }
  }

  if(!col->temperature_output_is_off())
  {
    if(col->getWriteMethod()=="pvtk")
    {
     if(!(col->getTemperatureOutputTag()).empty())
       temperaturewritebuffer=newArr4(float,grid->getNZN()-3, grid->getNYN()-3, grid->getNXN()-3,6);
    }
  }

  if (verbosity)
    dprintf("(6) initialized buffers for output to vtk files.");

  cq = SaveDirName + "/ConservedQuantities.txt";
  if (myrank == 0) 
  {
    ofstream my_file(cq.c_str());
    my_file.close();
  }

  if (verbosity)
    dprintf("(7) initialized ConservedQuantities.txt file.");

  rho = new double[ns];
  Ke = new double[ns];
  BulkEnergy = new double[ns];
  momentum = new double[ns];
  Qdel = new double[ns];
  Count = new double[ns];
  Qrep = new double[ns];
  Qexo = new double[ns];
  
  // put charges at zero to avoid random values
  for (int i = 0; i < ns; i++)
  {
    Qdel[i] = 0.;
    Count[i] = 0.;
    Qrep[i] = 0.;
    Qexo[i] = 0.;
  }

  //this initialization if for the 3D physical space 
  #ifdef USE_CATALYST_LEGACY
  Adaptor_legacy::Initialize(col, \
		  (int)(grid->getXstart()/grid->getDX()), \
		  (int)(grid->getYstart()/grid->getDY()), \
		  (int)(grid->getZstart()/grid->getDZ()), \
		  grid->getNXN(),
		  grid->getNYN(),
		  grid->getNZN(),
		  grid->getDX(),
		  grid->getDY(),
		  grid->getDZ());
  #endif

  my_clock = new Timing(myrank);

  if (verbosity)
    dprintf("(8) DONE.");

  return 0;
}


/*  ------------------- */
/*!  Calculate Moments  */
/*  ------------------- */
void c_Solver::CalculateMoments() {

  timeTasks_set_main_task(TimeTasks::MOMENTS);

  pad_particle_capacities();

  if (verbosity)
    dprintf("(1) First define the method to compute moments");

  if(Parameters::get_VECTORIZE_MOMENTS())
  {
    if (verbosity)
      dprintf("(2a) Particles are sorted by mesh cell = use vectorized moments");
    switch(Parameters::get_MOMENTS_TYPE())
    {
      case Parameters::SoA:
        convertParticlesToSoA();
        sortParticles();
        EMf->sumMoments_vectorized(part);
	if (verbosity)
          dprint("(2a) SoA method");
        break;
      case Parameters::AoS:
        convertParticlesToAoS();
        sortParticles();
        EMf->sumMoments_vectorized_AoS(part);
       	if (verbosity)
          dprint("(2a) Aos method");
        break;
      default:
        unsupported_value_error(Parameters::get_MOMENTS_TYPE());
    }
  }
  else
  {
    if (verbosity)
      dprintf("(2b) Particles are not sorted = use NON-vectorized moments");
    if(Parameters::get_SORTING_PARTICLES())
    {
      sortParticles();
      if (verbosity)
        dprintf("(2c) Particles are not sorted = sort particles by mesh cell");
    }
    switch(Parameters::get_MOMENTS_TYPE())
    {
      case Parameters::SoA:
        EMf->setZeroPrimaryMoments();
        convertParticlesToSoA();
        EMf->sumMoments(part);
	if (verbosity)
          dprint("(2b) SoA method");
        break;
      case Parameters::AoS:
        EMf->setZeroPrimaryMoments();
        convertParticlesToAoS();
        EMf->sumMoments_AoS(part);
	if (verbosity)
          dprint("(2b) AoS method");
        break;
      case Parameters::AoSintr:
        EMf->setZeroPrimaryMoments();
        convertParticlesToAoS();
        EMf->sumMoments_AoS_intr(part);
	if (verbosity)
          dprint("(2b) AoSintr method");
        break;
      default:
        unsupported_value_error(Parameters::get_MOMENTS_TYPE());
    }
  }

  EMf->setZeroDerivedMoments();
  EMf->sumOverSpecies();
  EMf->ConstantChargePlanet(col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
  
  if (verbosity) 
    dprintf("(3) Sum moments over species and apply constant charge inside planet");

  EMf->interpDensitiesN2C();
  
  if (verbosity) 
    dprintf("(4) Interpolate density from nodes to centers of mesh");
  
  EMf->calculateHatFunctions();

  if (verbosity) 
    dprintf("(5) Calculate hat functions (rho_hat, J_hat) for implicit solver");
  
}


/*  -------------- */
/*!  Calculate E  */
/*  -------------- */
void c_Solver::CalculateField(int cycle) {
  timeTasks_set_main_task(TimeTasks::FIELDS);
  EMf->calculateE(cycle);
  if (verbosity) 
    dprintf("(1) Calculate E-field with implicit method");
}


/*  -------------- */
/*!  Calculate B  */
/*  -------------- */
void c_Solver::CalculateB() {
  timeTasks_set_main_task(TimeTasks::FIELDS);
  EMf->calculateB();
  if (verbosity) 
    dprintf("(1) Calculate B-field from rot(E)");
}


/*  -------------- */
/*!  Particle mover */
/*  -------------- */
bool c_Solver::ParticlesMover(int cycle) {

  timeTasks_set_main_task(TimeTasks::PARTICLES);

  EMf->set_fieldForPcls();
  pad_particle_capacities();

  if (verbosity)
    dprintf("(1) Initialize fields for particles mover");

  bool applyCollisions = (col->getcollisionProcesses()) && (cycle % col->getcollStepSkip() == 0);

  for (int i = 0; i < ns; i++)
  {
    // should merely pass EMf->get_fieldForPcls() rather than EMf.
    // use the Predictor Corrector scheme to move particles
    switch(Parameters::get_MOVER_TYPE())
    {
      case Parameters::SoA:
        part[i].mover_PC(EMf);
        break;
      case Parameters::AoS:
        part[i].mover_PC_AoS(EMf);
        break;
      case Parameters::AoS_Relativistic:
        part[i].mover_PC_AoS_Relativistic(EMf);
        break;
      case Parameters::AoSintr:
        part[i].mover_PC_AoS_vec_intr(EMf);
        break;
      case Parameters::AoSvec:
        part[i].mover_PC_AoS_vec(EMf);
        break;
      default:
        unsupported_value_error(Parameters::get_MOVER_TYPE());
    }
 
    if (verbosity)
      dprintf("(2) PC pcls mover done. Species %d",i);
   
    if (applyCollisions)
    { 
      colls->Collide(i, part, col, EMf);
      if (verbosity)
        dprintf("(2a) Pcls collsion with neutrals. Species %d",i);
    }

    if ( (i>=col->getNs_sw()) and (col->getAddExosphere()) )
    {
      int i_pl = i-col->getNs_sw();
      Qexo[i] = part[i].AddIonizedExosphere(i_pl);
      if (verbosity and myrank==0){
        dprintf("(2b) Injection of particles due to photoionization. Species %d",i);
        dprintf("(2c) Ionized Exosph Injection the total Q(is=%d) is = %e",i,Qexo[i]);
      }
    }
  }
    
  if (applyCollisions)
  {
    colls->createIonizedParticles(part);
    if (verbosity)
      dprintf("(3a) Injection of particles due to electron impact ionization.");
  }

  double Qrm=0., Count_plus=0., Count_mins=0.;
    
  if (col->getNonTrivialBCPlanet()) 
  {
    for (int i=0; i < ns; i++)
    {
      Count[i] = part[i].rotateAndCountParticlesInsideSphere(cycle,col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
      if (Count[i]>0)  Count_plus += Count[i];
      if (Count[i]<0)  Count_mins += Count[i];
      if ( col->getVerbose() and Count[i]>0 and myrank==0) 
        dprintf("(4a) BC planet surface RotateAndCount-> the total Q(is=%d) counted is = %e",i,Count[i]);
    }
    Qrm = std::min(Count_plus,-Count_mins);
  }
  else 
    Qrm = (double) INT_MAX;

  for (int i=0; i < ns; i++) 
  {
    Qdel[i] = part[i].deleteParticlesInsideSphere(cycle,Qrm,col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
    if ( col->getVerbose() and fabs(Qdel[i])>0 and myrank==0)
      dprintf("(4b) BC planet surface: Delete-> the total Q(is=%d) removed is = %e",i,Qdel[i]);
  }

  for (int i=0; i < ns; i++)
  { 
    Qrep[i] = part[i].repopulate_particles(EMf);
    if (verbosity)
      dprintf("(5a) BC box walls (remove and inject): the total Q(is=%d) removed is = %e",i,Qrep[i]);

    // pcls communication - qui ci sta il problema su BCXright da capire...
    part[i].separate_and_send_particles();
    //part[i].communicate_particles();
    part[i].recommunicate_particles_until_done(1);
    if (verbosity)
      dprintf("(6a) Communication pcls different MPI procs DONE. Species %d",i);
  }

  return (false);

}


/*  -------------- */
/*!  Write Output */
/*  -------------- */
void c_Solver::WriteOutput(int cycle) {

  WriteConserved(cycle);

  #ifdef NO_HDF5
    eprintf("The selected output option must be compiled with HDF5");
  #else
    if (restart_cycle>0 && cycle%restart_cycle==0){
      convertParticlesToSynched();
      fetch_outputWrapperFPP().append_restart(cycle);
      if (verbosity) dprintf("(2) Dump all particles in restart*.hdf. Done.");
    }
  #endif  

  if (col->getWriteMethod() == "nbcvtk")
  {
    if(!col->field_output_is_off() && (cycle%(col->getFieldOutputCycle()) == 0 || cycle == first_cycle) )
    {
      if(!(col->getFieldOutputTag()).empty())
      {
        if(fieldreqcounter>0)
	{
	  //MPI_Waitall(fieldreqcounter,&fieldreqArr[0],&fieldstsArr[0]);\
	  // why above is commented ? ./job --> can merge the two if 
	  for(int si=0; si<fieldreqcounter; si++)
	  {
	    int error_code = MPI_File_write_all_end(fieldfhArr[si],&fieldwritebuffer[si][0][0][0],&fieldstsArr[si]);//fieldstsArr[si].MPI_ERROR;
	    if (error_code != MPI_SUCCESS) 
	    {
              char error_string[100];
              int length_of_error_string, error_class;
              MPI_Error_class(error_code, &error_class);
              MPI_Error_string(error_class, error_string, &length_of_error_string);
	      eprintf("MPI_Waitall error at field output cycle %d  %d  %s\n",cycle, si, error_string); // changed this from dprintf to printf ./job
            }
	    else
	    {
              MPI_File_close(&(fieldfhArr[si]));
            }
	  }
	}
	fieldreqcounter = WriteFieldsVTKNonblk(grid, EMf, col, vct,cycle,fieldwritebuffer,fieldreqArr,fieldfhArr);
      }
      if (verbosity) dprintf("(2a) nbcvtk done with the fields output.");

      #ifdef USE_CATALYST_LEGACY
        Adaptor_legacy::CoProcess(col->getDt()*cycle, cycle, EMf); //questo deve sta qui??? ./job
      #endif

      if(!(col->getMomentsOutputTag()).empty())
      {
        if(momentreqcounter>0)
	{
	//MPI_Waitall(momentreqcounter,&momentreqArr[0],&momentstsArr[0]);
	// same as above
	  for(int si=0;si< momentreqcounter;si++)
	  {
            int error_code = MPI_File_write_all_end(momentfhArr[si],&momentwritebuffer[si][0][0],&momentstsArr[si]);//momentstsArr[si].MPI_ERROR;
            if (error_code != MPI_SUCCESS) 
	    {
	      char error_string[100];
	      int length_of_error_string, error_class;
	      MPI_Error_class(error_code, &error_class);
	      MPI_Error_string(error_class, error_string, &length_of_error_string);
	      dprintf("MPI_Waitall error at moments output cycle %d  %d %s\n",cycle, si, error_string);
	    }
	    else
	    {
	      MPI_File_close(&(momentfhArr[si]));
	    }
	  }
	}
	momentreqcounter = WriteMomentsVTKNonblk(grid, EMf, col, vct,cycle,momentwritebuffer,momentreqArr,momentfhArr);
      }
      if (verbosity) dprintf("(2b) nbcvtk done with the moments output.");
    }
  }
  
  else if (col->getWriteMethod() == "pvtk")
  {
    if(!col->field_output_is_off() && (cycle%(col->getFieldOutputCycle()) == 0 || cycle == first_cycle) )
    {
      if(!(col->getFieldOutputTag()).empty())
        WriteFieldsVTK(grid, EMf, col, vct, col->getFieldOutputTag() ,cycle, fieldwritebuffer);//check this is E, B
      if (verbosity) dprintf("(2a) pvtk done with the fields output.");

      if(!(col->getMomentsOutputTag()).empty())
        WriteMomentsVTK(grid, EMf, col, vct, col->getMomentsOutputTag() ,cycle, momentwritebuffer);// and check this as Je0, Ji1 rhoe0 etc.
      if (verbosity) dprintf("(2b) pvtk done with the moments output.");
    }

    if(!col->spectra_output_is_off() && (cycle%(col->getSpectraOutputCycle()) == 0 || cycle == first_cycle) )
    {		  
      if(!(col->getSpectraOutputTag()).empty())
        WriteSpectraVTK(grid, part, EMf, col, vct, col->getSpectraOutputTag(),cycle, spectrawritebuffere, spectrawritebufferi);
      if (verbosity) dprintf("(2c) pvtk done with the spectra output.");
    }

    if(!col->temperature_output_is_off() && (cycle%(col->getTemperatureOutputCycle()) == 0 || cycle == first_cycle) )
    {
      if(!(col->getTemperatureOutputTag()).empty())
        WriteTemperatureVTK(grid, EMf, col, vct, col->getTemperatureOutputTag() ,cycle, temperaturewritebuffer);
      if (verbosity) dprintf("(2d) pvtk done with the Temperature output.");
    }
  }

  else if (col->getWriteMethod() == "phdf5")
  {
    if (!col->field_output_is_off() && cycle%(col->getFieldOutputCycle())==0 || cycle == first_cycle)
      WriteOutputParallel(grid, EMf, part, col, vct, cycle);
    if (verbosity) dprintf("(2a) phdf5 done with the fields and moments output (same file).");
  }

  else if (col->getWriteMethod() == "shdf5")
  {
    if (!col->field_output_is_off() && cycle%(col->getFieldOutputCycle())==0 || cycle == first_cycle)
    {
      if(!(col->getFieldOutputTag()).empty())
        fetch_outputWrapperFPP().append_output((col->getFieldOutputTag()).c_str(), cycle);
      if (verbosity) dprintf("(2a) shdf5 done with the fields output.");
      if(!(col->getMomentsOutputTag()).empty())
        fetch_outputWrapperFPP().append_output((col->getMomentsOutputTag()).c_str(), cycle);
      if (verbosity) dprintf("(2b) shdf5 done with the moments output.");
    }
   }

  else
  {
      warning_printf("Invalid output option. Options are: phdf5, shdf5, pvtk, nbcvtk");
      invalid_value_error(col->getWriteMethod().c_str());
  }

}


/*  -------------- */
/*!  Write Conserved */
/*  -------------- */
void c_Solver::WriteConserved(int cycle) 
{
  if(col->getDiagnosticsOutputCycle() > 0 && cycle % col->getDiagnosticsOutputCycle() == 0) 
  {
    if (verbosity) 
      dprintf("(1) START Write ConservedQuantities.txt");
    Eenergy = EMf->getEenergy();
    Benergy = EMf->getBenergy();
    TOTenergy = 0.0;
    TOTmomentum = 0.0;
    for (int is = 0; is < ns; is++) 
    {
      rho[is]= part[is].getRho();
      Ke[is] = part[is].getKe();
      BulkEnergy[is] = EMf->getBulkEnergy(is);
      TOTenergy += Ke[is];
      momentum[is] = part[is].getP();
      TOTmomentum += momentum[is];
    }
    if (myrank == (nprocs-1)) 
    {
      ofstream my_file(cq.c_str(), fstream::app);
      if(cycle==0)
      { 
        my_file << "#cycle" << "\t" << "Total_Energy" << "\t" << "Momentum" << "\t" << "Eenergy" << "\t" << "Benergy" << "\t" << "Kenergy" ;
	for (int is = 0; is < ns; is++) my_file << "Kenergy("<<is<<")" << "\t";
	for (int is = 0; is < ns; is++) my_file << "BulkEnergy("<<is<<")" << "\t";
	for (int is = 0; is < ns; is++) my_file << "Rho("<<is<<")" << "\t";
	for (int is = 0; is < ns; is++) my_file << "Qrep("<<is<<")" << "\t"; 
	for (int is = 0; is < ns; is++) my_file << "Qsphere("<<is<<")" << "\t"; 
	for (int is = 0; is < ns; is++) my_file << "Qexo("<<is<<")" << "\t"; 
	my_file << endl;
      }
      my_file << cycle << "\t" << (Eenergy + Benergy + TOTenergy) << "\t" << TOTmomentum << "\t" << Eenergy << "\t" << Benergy << "\t" << TOTenergy;
      for (int is = 0; is < ns; is++) my_file << "\t" << Ke[is];
      for (int is = 0; is < ns; is++) my_file << "\t" << BulkEnergy[is];
      for (int is = 0; is < ns; is++) my_file << "\t" << rho[is];
      for (int is = 0; is < ns; is++) my_file << "\t" << Qrep[is];
      for (int is = 0; is < ns; is++) my_file << "\t" << Qdel[is];
      for (int is = 0; is < ns; is++) my_file << "\t" << Qexo[is];
      my_file << endl;
      my_file.close();
    }
    if (verbosity) 
      dprintf("(2) DONE Write ConservedQuantities.txt");
  }
}

/*
void c_Solver::WriteVelocityDistribution(int cycle)
{
  // Velocity distribution
  //if(cycle % col->getVelocityDistributionOutputCycle() == 0)
  {
    for (int is = 0; is < ns; is++) {
      double maxVel = part[is].getMaxVelocity();
      long long *VelocityDist = part[is].getVelocityDistribution(nDistributionBins, maxVel);
      if (myrank == 0) {
        ofstream my_file(ds.c_str(), fstream::app);
        my_file << cycle << "\t" << is << "\t" << maxVel;
        for (int i = 0; i < nDistributionBins; i++)
          my_file << "\t" << VelocityDist[i];
        my_file << endl;
        my_file.close();
      }
      delete [] VelocityDist;
    }
  }
}

// This seems to record values at a grid of sample points
//
void c_Solver::WriteVirtualSatelliteTraces()
{
  if(ns <= 2) return;
  assert_eq(ns,4);

  ofstream my_file(cqsat.c_str(), fstream::app);
  const int nx0 = grid->get_nxc_r();
  const int ny0 = grid->get_nyc_r();
  const int nz0 = grid->get_nzc_r();
  for (int isat = 0; isat < nsat; isat++) {
    for (int jsat = 0; jsat < nsat; jsat++) {
      for (int ksat = 0; ksat < nsat; ksat++) {
        int index1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
        int index2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
        int index3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
        my_file << EMf->getBx(index1, index2, index3) << "\t" << EMf->getBy(index1, index2, index3) << "\t" << EMf->getBz(index1, index2, index3) << "\t";
        my_file << EMf->getEx(index1, index2, index3) << "\t" << EMf->getEy(index1, index2, index3) << "\t" << EMf->getEz(index1, index2, index3) << "\t";
        my_file << EMf->getJxs(index1, index2, index3, 0) + EMf->getJxs(index1, index2, index3, 2) << "\t" << EMf->getJys(index1, index2, index3, 0) + EMf->getJys(index1, index2, index3, 2) << "\t" << EMf->getJzs(index1, index2, index3, 0) + EMf->getJzs(index1, index2, index3, 2) << "\t";
        my_file << EMf->getJxs(index1, index2, index3, 1) + EMf->getJxs(index1, index2, index3, 3) << "\t" << EMf->getJys(index1, index2, index3, 1) + EMf->getJys(index1, index2, index3, 3) << "\t" << EMf->getJzs(index1, index2, index3, 1) + EMf->getJzs(index1, index2, index3, 3) << "\t";
        my_file << EMf->getRHOns(index1, index2, index3, 0) + EMf->getRHOns(index1, index2, index3, 2) << "\t";
        my_file << EMf->getRHOns(index1, index2, index3, 1) + EMf->getRHOns(index1, index2, index3, 3) << "\t";
      }}}
  my_file << endl;
  my_file.close();
}
*/

/*
void c_Solver::WriteParticles(int cycle)
{
#ifndef NO_HDF5
  if(col->particle_output_is_off() || cycle%(col->getParticlesOutputCycle())!=0) return;

  // this is a hack
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSynched();

  fetch_outputWrapperFPP().append_output((col->getPclOutputTag()).c_str(), cycle, 0);//"position + velocity + q "
#endif
}

void c_Solver::WriteTestParticles(int cycle)
{
#ifndef NO_HDF5
  if(nstestpart == 0 || col->testparticle_output_is_off() || cycle%(col->getTestParticlesOutputCycle())!=0) return;

  // this is a hack
  for (int i = 0; i < nstestpart; i++)
    testpart[i].convertParticlesToSynched();

  fetch_outputWrapperFPP().append_output("testpartpos + testpartvel+ testparttag", cycle, 0); // + testpartcharge
#endif
}
*/

// This needs to be separated into methods that save particles
// and methods that save field data
//
void c_Solver::Finalize() {
  if (verbosity) 
    dprintf("(1) START Finalize");
  if (col->getCallFinalize() && Parameters::get_doWriteOutput())
  {
    #ifndef NO_HDF5
    convertParticlesToSynched();
    fetch_outputWrapperFPP().append_restart((col->getNcycles() + first_cycle) - 1);
    #endif
  }
  my_clock->stopTiming();
  if (verbosity) 
    dprintf("(2) Stop time profiling. Finalize DONE.");
}

void c_Solver::sortParticles() {

  for(int species_idx=0; species_idx<ns; species_idx++)
    part[species_idx].sort_particles_serial();

}

void c_Solver::pad_particle_capacities()
{
  for (int i = 0; i < ns; i++)
    part[i].pad_capacities();

  for (int i = 0; i < nstestpart; i++)
    testpart[i].pad_capacities();
}

// convert particle to struct of arrays (assumed by I/O)
void c_Solver::convertParticlesToSoA()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSoA();
}

// convert particle to array of structs (used in computing)
void c_Solver::convertParticlesToAoS()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToAoS();
}

// convert particle to array of structs (used in computing)
void c_Solver::convertParticlesToSynched()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSynched();

  for (int i = 0; i < nstestpart; i++)
    testpart[i].convertParticlesToSynched();
}


int c_Solver::LastCycle() {
    return (col->getNcycles() + first_cycle);
}
