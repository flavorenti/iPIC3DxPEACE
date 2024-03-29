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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "input_array.h"
#include "ipichdf5.h"
#include "Collective.h"
#include "ConfigFile.h"
#include "limits.h" // for INT_MAX
#include "MPIdata.h"
#include "debug.h"
#include "asserts.h" // for assert_ge
#include "string.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>

// order must agree with Enum in Collective.h
static const char *enumNames[] =
{
  "default",
  "initial",
  "final",
  // used by ImplSusceptMode
  "explPredict",
  "implPredict",
  // marker for last enumerated symbol of this class
  "NUMBER_OF_ENUMS",
  "INVALID_ENUM"
};

int Collective::read_enum_parameter(const char* option_name, const char* default_value,
  const ConfigFile& config)
{
  string enum_name = config.read < string >(option_name,default_value);
  // search the list (could use std::map)
  //
  for(int i=0;i<NUMBER_OF_ENUMS;i++)
  {
    if(!strcmp(enum_name.c_str(),enumNames[i]))
      return i;
  }
  // could not find enum, so issue error and quit.
  if(!MPIdata::get_rank())
  {
    eprintf("in input file %s there is an invalid option %s\n",
      inputfile.c_str(), enum_name.c_str());
  }
  MPIdata::exit(1);
  // this is a better way
  return INVALID_ENUM;
}

const char* Collective::get_name_of_enum(int in)
{
  assert_ge(in, 0);
  assert_lt(in, NUMBER_OF_ENUMS);
  return enumNames[in];
}

/*! Read the input file from text file and put the data in a collective wrapper: if it's a restart read from input file basic sim data and load particles and EM field from restart file */
void Collective::ReadInput(string inputfile) {
  using namespace std;
  int test_verbose;
  // Loading the input file 
  ConfigFile config(inputfile);
  // the following variables are ALWAYS taken from inputfile, even if restarting 
  {
    /* first load SIMULATION PARAMETERS */	  
    SaveDirName       = config.read < string > ("SaveDirName","data");
    RestartDirName    = config.read < string > ("RestartDirName","data");
    Case              = config.read<string>("Case");
    wmethod           = config.read<string>("WriteMethod","pvtk");
    SimName           = config.read<string>("SimulationName");

    /* second load NUMERICAL PARAMETERS */	  
    dt                = config.read < double >("dt");
    ncycles           = config.read < int >("ncycles");
    th                = config.read < double >("th",1.0);
    c                 = config.read < double >("c",1.0);
    Smooth            = config.read < double >("Smooth",1.0);
    SmoothNiter       = config.read < int >("SmoothNiter",6);
    PoissonCorrection = config.read<string>("PoissonCorrection","yes");
    PoissonCorrectionCycle = config.read<int>("PoissonCorrectionCycle",10);
    Lx                = config.read < double >("Lx",10.0);
    Ly                = config.read < double >("Ly",10.0);
    Lz                = config.read < double >("Lz",10.0);
    nxc               = config.read < int >("nxc",64);
    nyc               = config.read < int >("nyc",64);
    nzc               = config.read < int >("nzc",64);
    XLEN              = config.read < int >("XLEN",1);
    YLEN              = config.read < int >("YLEN",1);
    ZLEN              = config.read < int >("ZLEN",1);
    PERIODICX         = config.read < bool >("PERIODICX",true);
    PERIODICY         = config.read < bool >("PERIODICY",true);
    PERIODICZ         = config.read < bool >("PERIODICZ",true);
    PERIODICX_P       = config.read < bool >("PERIODICX_P",PERIODICX);
    PERIODICY_P       = config.read < bool >("PERIODICY_P",PERIODICY);
    PERIODICZ_P       = config.read < bool >("PERIODICZ_P",PERIODICZ);
    ns                = config.read < int >("ns",2);
    nstestpart        = config.read < int >("nsTestPart", 0); //TBR?
    NpMaxNpRatio      = config.read < double >("NpMaxNpRatio",1.5); //TBR?
    assert_ge(NpMaxNpRatio, 1.);
    verbose           = config.read < bool > ("verbose",true);
    delta             = config.read < double >("delta",0.5); //TBR
    CGtol             = config.read < double >("CGtol",1e-3);
    GMREStol          = config.read < double >("GMREStol",1e-3);
    NiterMover        = config.read < int >("NiterMover",8);
    int ns_tot =ns+nstestpart;
    npcelx            = new int[ns_tot];
    npcely            = new int[ns_tot];
    npcelz            = new int[ns_tot];
    array_int npcelx0 = config.read < array_int > ("npcelx");
    array_int npcely0 = config.read < array_int > ("npcely");
    array_int npcelz0 = config.read < array_int > ("npcelz");
    npcelx[0] = npcelx0.a;
    npcely[0] = npcely0.a;
    npcelz[0] = npcelz0.a;
    if (ns_tot > 1) {
      npcelx[1] = npcelx0.b;
      npcely[1] = npcely0.b;
      npcelz[1] = npcelz0.b;
    }
    if (ns_tot > 2) {
      npcelx[2] = npcelx0.c;
      npcely[2] = npcely0.c;
      npcelz[2] = npcelz0.c;
    }
    if (ns_tot > 3) {
      npcelx[3] = npcelx0.d;
      npcely[3] = npcely0.d;
      npcelz[3] = npcelz0.d;
    }
    if (ns_tot > 4) {
      npcelx[4] = npcelx0.e;
      npcely[4] = npcely0.e;
      npcelz[4] = npcelz0.e;
    }
    if (ns_tot > 5) {
      npcelx[5] = npcelx0.f;
      npcely[5] = npcely0.f;
      npcelz[5] = npcelz0.f;
    }
    if (ns_tot > 6) {
      npcelx[6] = npcelx0.g;
      npcely[6] = npcely0.g;
      npcelz[6] = npcelz0.g;
    }
    if (ns_tot > 7) {
      npcelx[7] = npcelx0.h;
      npcely[7] = npcely0.h;
      npcelz[7] = npcelz0.h;
    }
    if (ns_tot > 8) {
      npcelx[8] = npcelx0.i;
      npcely[8] = npcely0.i;
      npcelz[8] = npcelz0.i;
    }
    if (ns_tot > 9) {
      npcelx[9] = npcelx0.j;
      npcely[9] = npcely0.j;
      npcelz[9] = npcelz0.j;
    }
    if (ns_tot > 10) {
      npcelx[10] = npcelx0.k;
      npcely[10] = npcely0.k;
      npcelz[10] = npcelz0.k;
    }
    if (ns_tot > 11) {
      npcelx[11] = npcelx0.l;
      npcely[11] = npcely0.l;
      npcelz[11] = npcelz0.l;
    }

    /* third load SOLAR WIND PARAMETERS */	  
    B0x               = config.read <double>("B0x",0.0);
    B0y               = config.read <double>("B0y",0.0);
    B0z               = config.read <double>("B0z",0.0);
    ns_sw             = config.read < int >("ns_sw",2);
    rhoINIT           = new double[ns_sw];
    array_double rhoINIT0 = config.read < array_double > ("rhoINIT");
    rhoINIT[0] = rhoINIT0.a;
    if (ns_sw > 1)
      rhoINIT[1] = rhoINIT0.b;
    if (ns_sw > 2)
      rhoINIT[2] = rhoINIT0.c;
    if (ns_sw > 3)
      rhoINIT[3] = rhoINIT0.d;
    if (ns_sw > 4)
      rhoINIT[4] = rhoINIT0.e;
    if (ns_sw > 5)
      rhoINIT[5] = rhoINIT0.f;
    Vinj              = config.read < double >("Vinj",0.0);  //TBR

    /* fourth load PLANET PARAMETERS */	  
    B1z               = config.read <double>("B1z",0.0);
    x_center          = config.read < double >("x_center",5.0);
    y_center          = config.read < double >("y_center",5.0);
    z_center          = config.read < double >("z_center",5.0);
    L_square          = config.read < double >("L_square",5.0);
    PlanetOffset      = config.read <double>("PlanetOffset",0.0);
    ns_pl             = config.read < int >("ns_pl",0);
    AddExosphere      = config.read < int >("AddExosphere",0);
    Rmax              = config.read < double >("Rmax",3.);
    Wfact = new double[ns_pl];
    nSurf = new double[ns_pl];
    hExo = new double[ns_pl];
    fExo = new double[ns_pl];
    array_double Wfact0;
    array_double nSurf0;
    array_double hExo0;
    array_double fExo0;
    if (ns_pl > 0)
    {
      Wfact0 = config.read < array_double > ("Weight_factor");
      nSurf0 = config.read < array_double > ("nSurf");
      hExo0 = config.read < array_double > ("hExo");
      fExo0 = config.read < array_double > ("fExo");
      Wfact[0] = Wfact0.a;
      nSurf[0] = nSurf0.a;
      hExo[0] = hExo0.a;
      fExo[0] = fExo0.a;
    }
    if (ns_pl > 1)
    {
      Wfact[1] = Wfact0.b;
      nSurf[1] = nSurf0.b;
      hExo[1]  = hExo0.b;
      fExo[1] = fExo0.b;
    }
    if (ns_pl > 2)
    {
      Wfact[2] = Wfact0.c;
      nSurf[2] = nSurf0.c;
      hExo[2]  = hExo0.c;
      fExo[2] = fExo0.c;
    }
    if (ns_pl > 3)
    {
      Wfact[3] = Wfact0.d;
      nSurf[3] = nSurf0.d;
      hExo[3]  = hExo0.d;
      fExo[3] = fExo0.d;
    }
    if (ns_pl > 4)
    {
      Wfact[4] = Wfact0.e;
      nSurf[4] = nSurf0.e;
      hExo[4]  = hExo0.e;
      fExo[4] = fExo0.e;
    }
    if (ns_pl > 5)
    {
      Wfact[5] = Wfact0.f;
      nSurf[5] = nSurf0.f;
      hExo[5]  = hExo0.f;
      fExo[5] = fExo0.f;
    }
    collisionProcesses = config.read < bool >("collisionProcesses", 0);
    xSec               = config.read < double >("xSec", 8.82e-10);
    nCollProcesses     = config.read < int >("nCollProcesses", 3);
    nIoniColls         = config.read < int >("nIoniColls", 1);
    iSecElec           = config.read < int >("iSecElec", 2);
    iSecIon            = config.read < int >("iSecIon", 3);
    collStepSkip       = config.read< int >("collStepSkip", 1);
    E_th_el = new double[nCollProcesses];
    array_double E_th_el0;
    if ((ns_pl>0) and (nCollProcesses>0)) 
    { 
      E_th_el0 = config.read < array_double > ("E_th_el");
      E_th_el[0] = E_th_el0.a;
    }
    if (nCollProcesses > 1)
      E_th_el[1] = E_th_el0.b;
    if (nCollProcesses > 2)
      E_th_el[2] = E_th_el0.c;
    if (nCollProcesses > 3)
      E_th_el[3] = E_th_el0.d;
    if (nCollProcesses > 4)
      E_th_el[4] = E_th_el0.e;
    if (nCollProcesses > 5)
      E_th_el[5] = E_th_el0.f;

    /* fifth load the BOUNDARY PARAMETERS */
    yes_sal = config.read < int >("yes_sal",1);
    n_layers_sal = config.read < int >("n_layers_sal",3);
    NewPclInit = config.read < int >("NewPclInit",1); 
    NonTrivialBCPlanet = config.read < int >("NonTrivialBCPlanet",1);
    bcPHIfaceXright = config.read < int >("bcPHIfaceXright",1);
    bcPHIfaceXleft  = config.read < int >("bcPHIfaceXleft",1);
    bcPHIfaceYright = config.read < int >("bcPHIfaceYright",1);
    bcPHIfaceYleft  = config.read < int >("bcPHIfaceYleft",1);
    bcPHIfaceZright = config.read < int >("bcPHIfaceZright",1);
    bcPHIfaceZleft  = config.read < int >("bcPHIfaceZleft",1);
    bcEMfaceXright = config.read < int >("bcEMfaceXright",2);
    bcEMfaceXleft  = config.read < int >("bcEMfaceXleft",2);
    bcEMfaceYright = config.read < int >("bcEMfaceYright",2);
    bcEMfaceYleft  = config.read < int >("bcEMfaceYleft",2);
    bcEMfaceZright = config.read < int >("bcEMfaceZright",2);
    bcEMfaceZleft  = config.read < int >("bcEMfaceZleft",2);
    bcPfaceXright = config.read < int >("bcPfaceXright",2);
    bcPfaceXleft  = config.read < int >("bcPfaceXleft",2);
    bcPfaceYright = config.read < int >("bcPfaceYright",2);
    bcPfaceYleft  = config.read < int >("bcPfaceYleft",2);
    bcPfaceZright = config.read < int >("bcPfaceZright",2);
    bcPfaceZleft  = config.read < int >("bcPfaceZleft",2);

    /* sixth load the OUTPUT PARAMETERS */
    FieldOutputCycle   = config.read < int >("FieldOutputCycle",100);
    FieldOutputTag     = config.read <string>("FieldOutputTag","");
    MomentsOutputTag   = config.read <string>("MomentsOutputTag","");
    TemperatureOutputCycle = config.read < int >("TemperatureOutputCycle",0);     
    TemperatureOutputTag   = config.read <string>("TemperatureOutputTag","");
    SpectraOutputCycle = config.read < int >("SpectraOutputCycle",0);
    SpectraOutputTag   = config.read <string>("SpectraOutputTag","");
    ParticlesOutputCycle = config.read < int >("ParticlesOutputCycle",0);
    ParticlesOutputTag =   config.read <string>("ParticlesOutputTag","");
    testPartFlushCycle = config.read < int >("TestParticlesOutputCycle",0); //TBR ?
    RestartOutputCycle = config.read < int >("RestartOutputCycle",5000);
    RemoveParticlesOutputCycle = config.read < int >("RemoveParticlesOutputCycle",0);
    DiagnosticsOutputCycle = config.read < int >("DiagnosticsOutputCycle",10);
    DeltaX = config.read < int >("DeltaX",1);
    DeltaY = config.read < int >("DeltaY",1);
    DeltaZ = config.read < int >("DeltaZ",1);
    Estarti = config.read < double >("Estarti",-1);
    Eendi = config.read < double >("Eendi",1);
    dEi = config.read < double >("dEi",0.5);
    Estarte = config.read < double >("Estarte",-1);
    Eende = config.read < double >("Eende",1);
    dEe = config.read < double >("dEe",0.5);
    CallFinalize = config.read < bool >("CallFinalize", true); //deprecated TBR
    ParaviewScriptPath =config.read <string>("ParaviewScriptPath", "");
  }

  last_cycle = -1;

  // now read all the vth, v0 in a good way, merging sw and planet data
  ns_all  = ns_sw+ns_pl;

  if (ns_all!=ns)
    eprintf("ERROR number of species: ns_sw + ns_pl different from ns");

  uth = new double[ns_all];
  vth = new double[ns_all];
  wth = new double[ns_all];
  u0 = new double[ns_all];
  v0 = new double[ns_all];
  w0 = new double[ns_all];
  qom = new double[ns_all];

  array_double uth0_sw = config.read < array_double > ("uth_sw");
  array_double vth0_sw = config.read < array_double > ("vth_sw");
  array_double wth0_sw = config.read < array_double > ("wth_sw");
  array_double u00_sw = config.read < array_double > ("u0_sw");
  array_double v00_sw = config.read < array_double > ("v0_sw");
  array_double w00_sw = config.read < array_double > ("w0_sw");
  array_double qom_sw = config.read < array_double > ("qom_sw");
  
  array_double uth0_pl;
  array_double vth0_pl;
  array_double wth0_pl;
  array_double u00_pl;
  array_double v00_pl;
  array_double w00_pl;
  array_double qom_pl;

  if (ns_pl>0)
  {
    uth0_pl = config.read < array_double > ("uth_pl");
    vth0_pl = config.read < array_double > ("vth_pl");
    wth0_pl = config.read < array_double > ("wth_pl");
    u00_pl = config.read < array_double > ("u0_pl");
    v00_pl = config.read < array_double > ("v0_pl");
    w00_pl = config.read < array_double > ("w0_pl");
    qom_pl = config.read < array_double > ("qom_pl");
  }

  uth[0] = uth0_sw.a;
  vth[0] = vth0_sw.a;
  wth[0] = wth0_sw.a;
  u0[0] = u00_sw.a;
  v0[0] = v00_sw.a;
  w0[0] = w00_sw.a;
  qom[0] = qom_sw.a;

  uth[1] = uth0_sw.b;
  vth[1] = vth0_sw.b;
  wth[1] = wth0_sw.b;
  u0[1] = u00_sw.b;
  v0[1] = v00_sw.b;
  w0[1] = w00_sw.b;
  qom[1] = qom_sw.b;

  if (ns_sw > 2) {
    uth[2] = uth0_sw.c;
    vth[2] = vth0_sw.c;
    wth[2] = wth0_sw.c;
    u0[2] = u00_sw.c;
    v0[2] = v00_sw.c;
    w0[2] = w00_sw.c;
    qom[2] = qom_sw.c;
  }
  if (ns_sw > 3) {
    uth[3] = uth0_sw.d;
    vth[3] = vth0_sw.d;
    wth[3] = wth0_sw.d;
    u0[3] = u00_sw.d;
    v0[3] = v00_sw.d;
    w0[3] = w00_sw.d;
    qom[3] = qom_sw.d;
  }
  if (ns_sw > 4) {
    uth[4] = uth0_sw.e;
    vth[4] = vth0_sw.e;
    wth[4] = wth0_sw.e;
    u0[4] = u00_sw.e;
    v0[4] = v00_sw.e;
    w0[4] = w00_sw.e;
    qom[4] = qom_sw.e;
  }
  if (ns_sw > 5) {
    uth[5] = uth0_sw.f;
    vth[5] = vth0_sw.f;
    wth[5] = wth0_sw.f;
    u0[5] = u00_sw.f;
    v0[5] = v00_sw.f;
    w0[5] = w00_sw.f;
    qom[5] = qom_sw.f;
  }
  if (ns_pl > 0) {
    uth[ns_sw] = uth0_pl.a;
    vth[ns_sw] = vth0_pl.a;
    wth[ns_sw] = wth0_pl.a;
    u0[ns_sw] = u00_pl.a;
    v0[ns_sw] = v00_pl.a;
    w0[ns_sw] = w00_pl.a;
    qom[ns_sw] = qom_pl.a;
  }
  if (ns_pl > 1) {
    uth[ns_sw+1] = uth0_pl.b;
    vth[ns_sw+1] = vth0_pl.b;
    wth[ns_sw+1] = wth0_pl.b;
    u0[ns_sw+1] = u00_pl.b;
    v0[ns_sw+1] = v00_pl.b;
    w0[ns_sw+1] = w00_pl.b;
    qom[ns_sw+1] = qom_pl.b;
  }
  if (ns_pl > 2) {
    uth[ns_sw+2] = uth0_pl.c;
    vth[ns_sw+2] = vth0_pl.c;
    wth[ns_sw+2] = wth0_pl.c;
    u0[ns_sw+2] = u00_pl.c;
    v0[ns_sw+2] = v00_pl.c;
    w0[ns_sw+2] = w00_pl.c;
    qom[ns_sw+2] = qom_pl.c;
  }
  if (ns_pl > 3) {
    uth[ns_sw+3] = uth0_pl.d;
    vth[ns_sw+3] = vth0_pl.d;
    wth[ns_sw+3] = wth0_pl.d;
    u0[ns_sw+3] = u00_pl.d;
    v0[ns_sw+3] = v00_pl.d;
    w0[ns_sw+3] = w00_pl.d;
    qom[ns_sw+3] = qom_pl.d;
  }
  if (ns_pl > 4) {
    uth[ns_sw+4] = uth0_pl.e;
    vth[ns_sw+4] = vth0_pl.e;
    wth[ns_sw+4] = wth0_pl.e;
    u0[ns_sw+4] = u00_pl.e;
    v0[ns_sw+4] = v00_pl.e;
    w0[ns_sw+4] = w00_pl.e;
    qom[ns_sw+4] = qom_pl.e;
  }
  if (ns_pl > 5) {
    uth[ns_sw+5] = uth0_pl.f;
    vth[ns_sw+5] = vth0_pl.f;
    wth[ns_sw+5] = wth0_pl.f;
    u0[ns_sw+5] = u00_pl.f;
    v0[ns_sw+5] = v00_pl.f;
    w0[ns_sw+5] = w00_pl.f;
    qom[ns_sw+5] = qom_pl.f;
  }

  // electric field in SW
  E0x = w0[0]*B0y - v0[0]*B0z;
  E0y =-u0[0]*B0z + w0[0]*B0x;
  E0z = u0[0]*B0y - v0[0]*B0x;


  // translates bcEM into BC for EM fields E and B
  bcEx[0] = bcEMfaceXright == 0 ? 2 : 1;   bcBx[0] = bcEMfaceXright == 0 ? 1 : 2;
  bcEx[1] = bcEMfaceXleft  == 0 ? 2 : 1;   bcBx[1] = bcEMfaceXleft  == 0 ? 1 : 2;
  bcEx[2] = bcEMfaceYright == 0 ? 1 : 2;   bcBx[2] = bcEMfaceYright == 0 ? 2 : 1;
  bcEx[3] = bcEMfaceYleft  == 0 ? 1 : 2;   bcBx[3] = bcEMfaceYleft  == 0 ? 2 : 1;
  bcEx[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBx[4] = bcEMfaceZright == 0 ? 2 : 1;
  bcEx[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBx[5] = bcEMfaceZleft  == 0 ? 2 : 1;

  bcEy[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBy[0] = bcEMfaceXright == 0 ? 2 : 1;
  bcEy[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBy[1] = bcEMfaceXleft  == 0 ? 2 : 1;
  bcEy[2] = bcEMfaceYright == 0 ? 2 : 1;   bcBy[2] = bcEMfaceYright == 0 ? 1 : 2;
  bcEy[3] = bcEMfaceYleft  == 0 ? 2 : 1;   bcBy[3] = bcEMfaceYleft  == 0 ? 1 : 2;
  bcEy[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBy[4] = bcEMfaceZright == 0 ? 2 : 1;
  bcEy[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBy[5] = bcEMfaceZleft  == 0 ? 2 : 1;

  bcEz[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBz[0] = bcEMfaceXright == 0 ? 2 : 1;
  bcEz[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBz[1] = bcEMfaceXleft  == 0 ? 2 : 1;
  bcEz[2] = bcEMfaceYright == 0 ? 1 : 1;   bcBz[2] = bcEMfaceYright == 0 ? 2 : 1;
  bcEz[3] = bcEMfaceYleft  == 0 ? 1 : 1;   bcBz[3] = bcEMfaceYleft  == 0 ? 2 : 1;
  bcEz[4] = bcEMfaceZright == 0 ? 2 : 1;   bcBz[4] = bcEMfaceZright == 0 ? 1 : 2;
  bcEz[5] = bcEMfaceZleft  == 0 ? 2 : 1;   bcBz[5] = bcEMfaceZleft  == 0 ? 1 : 2;

  #ifndef NO_HDF5 
  if (RESTART1) {               // you are restarting
    RestartDirName = config.read < string > ("RestartDirName","data");
    // ReadRestart(RestartDirName);
    restart_status = 1;

    hid_t file_id = H5Fopen((RestartDirName + "/restart0.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "couldn't open file: " << (RestartDirName + "/restart0.hdf").c_str() << endl;
      return;
    }

    hid_t dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);  // HDF 1.8.8
    herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
  }
#endif
}

bool Collective::field_output_is_off()const
{
  return (FieldOutputCycle <= 0);
}

bool Collective::spectra_output_is_off()const
{
  return (SpectraOutputCycle <= 0);
}

bool Collective::temperature_output_is_off()const
{
  return (TemperatureOutputCycle <= 0);
}

bool Collective::particle_output_is_off()const
{
  return getParticlesOutputCycle() <= 0;
}
bool Collective::testparticle_output_is_off()const
{
  return getTestParticlesOutputCycle() <= 0;
}


/*
*! Read the collective information from the RESTART file in HDF5 format
 * There are three restart status: restart_status = 0 ---> new inputfile
 * restart_status = 1 ---> RESTART and restart and result directories does not coincide
 * restart_status = 2 ---> RESTART and restart and result directories coincide *
int Collective::ReadRestart(string inputfile) { 

	*Routine deprecated F.Lavorenti*

#ifdef NO_HDF5
  eprintf("restart requires compiling with HDF5");
#else
  restart_status = 1;
  // hdf stuff 
  hid_t file_id;
  hid_t dataset_id;
  herr_t status;
  printf("\n Commencing restart with H5F commands. In collective.cpp. \n");
  // Open the setting file for the restart.
  file_id = H5Fopen((inputfile + "/settings.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: " << (inputfile + "/settings.hdf").c_str() << endl;
    return -1;
  }

  // read c
  dataset_id = H5Dopen2(file_id, "/collective/c", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c);
  status = H5Dclose(dataset_id);

  // read Lx 
  dataset_id = H5Dopen2(file_id, "/collective/Lx", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
  status = H5Dclose(dataset_id);
  // read Ly 
  dataset_id = H5Dopen2(file_id, "/collective/Ly", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
  status = H5Dclose(dataset_id);
  // read Lz 
  dataset_id = H5Dopen2(file_id, "/collective/Lz", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lz);
  status = H5Dclose(dataset_id);
  // read x_center
  dataset_id = H5Dopen2(file_id, "/collective/x_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x_center);
  status = H5Dclose(dataset_id);
  // read y_center
  dataset_id = H5Dopen2(file_id, "/collective/y_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y_center);
  status = H5Dclose(dataset_id);
  // read z_center
  dataset_id = H5Dopen2(file_id, "/collective/z_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z_center);
  status = H5Dclose(dataset_id);
  // read L_square
  dataset_id = H5Dopen2(file_id, "/collective/L_square", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L_square);
  status = H5Dclose(dataset_id);
  // read nxc
  dataset_id = H5Dopen2(file_id, "/collective/Nxc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxc);
  status = H5Dclose(dataset_id);
  // read nyc 
  dataset_id = H5Dopen2(file_id, "/collective/Nyc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyc);
  status = H5Dclose(dataset_id);
  // read nyc 
  dataset_id = H5Dopen2(file_id, "/collective/Nzc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzc);
  status = H5Dclose(dataset_id);
  // read ns
  dataset_id = H5Dopen2(file_id, "/collective/Ns", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ns);
  status = H5Dclose(dataset_id);
  //read number of test particles species
  dataset_id = H5Dopen2(file_id, "/collective/NsTestPart", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nstestpart);
  status = H5Dclose(dataset_id);

  *! Boundary condition information *
  // read EMfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceXleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceXleft);
  status = H5Dclose(dataset_id);
  // read EMfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceXright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceXright);
  status = H5Dclose(dataset_id);
  // read EMfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceYleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceYleft);
  status = H5Dclose(dataset_id);
  // read EMfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceYright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceYright);
  status = H5Dclose(dataset_id);
  // read EMfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceZleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceZleft);
  status = H5Dclose(dataset_id);
  // read EMfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceZright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceZright);
  status = H5Dclose(dataset_id);

  // read PHIfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceXleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceXleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceXright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceXright);
  status = H5Dclose(dataset_id);
  // read PHIfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceYleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceYleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceYright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceYright);
  status = H5Dclose(dataset_id);
  // read PHIfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceZleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceZleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceZright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceZright);
  status = H5Dclose(dataset_id);

  // read PfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceXleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceXleft);
  status = H5Dclose(dataset_id);
  // read PfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceXright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceXright);
  status = H5Dclose(dataset_id);
  // read PfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceYleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceYleft);
  status = H5Dclose(dataset_id);
  // read PfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceYright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceYright);
  status = H5Dclose(dataset_id);
  // read PfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceZleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceZleft);
  status = H5Dclose(dataset_id);
  // read PfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceZright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceZright);
  status = H5Dclose(dataset_id);
  // allocate fields depending on species
  npcelx = new int[ns+nstestpart];
  npcely = new int[ns+nstestpart];
  npcelz = new int[ns+nstestpart];
  qom = new double[ns+nstestpart];
  uth = new double[ns];
  vth = new double[ns];
  wth = new double[ns];
  u0 = new double[ns];
  v0 = new double[ns];
  w0 = new double[ns];
  // read data from species0, species 1, species2,...
  string *name_species = new string[ns];
  stringstream *ss = new stringstream[ns];
  string *name_testspecies;
  stringstream *testss;

  for (int i = 0; i < ns; i++) {
    ss[i] << i;
    name_species[i] = "/collective/species_" + ss[i].str() + "/";
  }
  if(nstestpart>0){
	  name_testspecies = new string[nstestpart];
	  testss = new stringstream[nstestpart];
	  for (int i = 0; i < nstestpart; i++) {
		  testss[i] << (i+ns);
		  name_testspecies[i] = "/collective/testspecies_" + testss[i].str() + "/";
	  }

	  pitch_angle = new double[nstestpart];
	  energy      = new double[nstestpart];
	  for (int i = 0; i < nstestpart; i++) {
	    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "pitch_angle").c_str(), H5P_DEFAULT); // HDF 1.8.8
	    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pitch_angle[i]);
	    status = H5Dclose(dataset_id);

	    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "energy").c_str(), H5P_DEFAULT); // HDF 1.8.8
	    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &energy[i]);
	    status = H5Dclose(dataset_id);
	  }
  }

  // npcelx for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcelx").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelx[i]);
    status = H5Dclose(dataset_id);
  }
  // npcely for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcely").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcely[i]);
    status = H5Dclose(dataset_id);
  }
  // npcelz for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcelz").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelz[i]);
    status = H5Dclose(dataset_id);
  }
  // qom for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "qom").c_str(), H5P_DEFAULT);  // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &qom[i]);
    status = H5Dclose(dataset_id);
  }

  //Test Particle
  // npcelx for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "Npcelx").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelx[i]);
    status = H5Dclose(dataset_id);
  }
  // npcely for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "Npcely").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcely[i]);
    status = H5Dclose(dataset_id);
  }
  // npcelz for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "Npcelz").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelz[i]);
    status = H5Dclose(dataset_id);
  }
  // qom for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "qom").c_str(), H5P_DEFAULT);  // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &qom[i]);
    status = H5Dclose(dataset_id);
  }


  *! not needed for restart * *
  for (int i = 0; i < ns; i++)
    uth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    vth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    wth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    u0[i] = 0.0;
  for (int i = 0; i < ns; i++)
    v0[i] = 0.0;
  for (int i = 0; i < ns; i++)
    w0[i] = 0.0;
  // verbose on
  //verbose = 1;

  // read collisionProcesses
  dataset_id = H5Dopen2(file_id, "/collective/colls/collisionProcesses", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &collisionProcesses);
  status = H5Dclose(dataset_id);
  // read xSec
  dataset_id = H5Dopen2(file_id, "/collective/colls/xSec", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &xSec);
  status = H5Dclose(dataset_id);
  // read iSecElec
  dataset_id = H5Dopen2(file_id, "/collective/colls/iSecElec", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &iSecElec);
  status = H5Dclose(dataset_id);
  // read iSecIon
  dataset_id = H5Dopen2(file_id, "/collective/colls/iSecIon", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &iSecIon);
  status = H5Dclose(dataset_id);
  // read nCollProcesses
  dataset_id = H5Dopen2(file_id, "/collective/colls/nCollProcesses", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nCollProcesses);
  status = H5Dclose(dataset_id);
  // allocate fields depending on nCollProcesses
  E_th_el = new double[nCollProcesses];
  // read data from E_th_el0, E_th_el1, E_th_el2,...
  string *nameEth = new string[nCollProcesses];
  stringstream *sColl = new stringstream[nCollProcesses];
  // read E_th_el
  for (int i = 0; i < ns; i++) {
    	sColl[i] << i;
    	nameEth[i] = "/collective/colls/E_th_el"+sColl[i].str();
	if (collisionProcesses){
    		dataset_id = H5Dopen2(file_id, (nameEth[i]).c_str(), H5P_DEFAULT); // HDF 1.8.8
    		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &E_th_el[i]);
    		status = H5Dclose(dataset_id);
  	}
  }
  // read nIoniColls
  dataset_id = H5Dopen2(file_id, "/collective/colls/nIoniColls", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nIoniColls);
  status = H5Dclose(dataset_id);
  // read collStepSkip
  dataset_id = H5Dopen2(file_id, "/collective/colls/collStepSkip", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &collStepSkip);
  status = H5Dclose(dataset_id);
  // read nNeutSpecies
  dataset_id = H5Dopen2(file_id, "/collective/colls/nNeutSpecies", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nNeutSpecies);
  status = H5Dclose(dataset_id);

  // allocate fields depending on nNeutSpecies
  nSurf = new double[nNeutSpecies];
  hExo = new double[nNeutSpecies];
  fExo = new double[nNeutSpecies];

  stringstream *sNeut = new stringstream[nNeutSpecies];
  for (int i = 0; i < nNeutSpecies; i++){
    sNeut[i] << i;
    // read nSurf
    dataset_id = H5Dopen2(file_id, ("/collective/colls/nSurf" + sNeut[i].str()).c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nSurf[i]);
    status = H5Dclose(dataset_id);
    // read hExo
    dataset_id = H5Dopen2(file_id, ("/collective/colls/hExo" + sNeut[i].str()).c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &hExo[i]);
    status = H5Dclose(dataset_id);
    // read fExo
    dataset_id = H5Dopen2(file_id, ("/collective/colls/fExo" + sNeut[i].str()).c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fExo[i]);
    status = H5Dclose(dataset_id);
  }


  // if RestartDirName == SaveDirName overwrite dt,Th,Smooth (append to old hdf files)
  if (RestartDirName == SaveDirName) {
    restart_status = 2;
    // read dt
    dataset_id = H5Dopen2(file_id, "/collective/Dt", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dt);
    status = H5Dclose(dataset_id);
    // read th 
    dataset_id = H5Dopen2(file_id, "/collective/Th", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &th);
    status = H5Dclose(dataset_id);
    // read Smooth
    dataset_id = H5Dopen2(file_id, "/collective/Smooth", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Smooth);
    status = H5Dclose(dataset_id);
    dataset_id = H5Dopen2(file_id, "/collective/SmoothNiter", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &SmoothNiter);
    status = H5Dclose(dataset_id);
  }

  status = H5Fclose(file_id);


  // read last cycle (not from settings, but from restart0.hdf)

  file_id = H5Fopen((inputfile + "/restart0.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: " << inputfile << endl;
    return -1;
  }

  dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  // deallocate
  delete[]name_species;
  delete[]ss;
  delete[]sColl;
  delete[]sNeut;
  if (collisionProcesses)
	  delete[]nameEth;
#endif
  return (0);
}
*/


void Collective::read_field_restart(
    const VCtopology3D* vct,
    const Grid* grid,
    arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
    arr3_double Ex, arr3_double Ey, arr3_double Ez,
    array4_double* rhons_, int ns)const
{
#ifdef NO_HDF5
  eprintf("Require HDF5 to read from restart file.");
#else
    const int nxn = grid->getNXN();
    const int nyn = grid->getNYN();
    const int nzn = grid->getNZN();
    if (vct->getCartesian_rank() == 0)
    {
      printf("LOADING EM FIELD FROM RESTART FILE in %s/restart.hdf\n",getRestartDirName().c_str());
    }

    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = getRestartDirName() + "/restart" + ss.str() + ".hdf";

    // hdf stuff
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[3];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      eprintf("Failed to open file: %s\n ", name_file.c_str());
    }

    //find the last cycle
    int lastcycle=0;
    dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &lastcycle);
    status = H5Dclose(dataset_id);

    // Bxn
    ss.str("");ss << "/fields/Bx/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    double *temp_storage = new double[dims_out[0] * dims_out[1] * dims_out[2]];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bxn[i][j][jj] = temp_storage[k++];


    status = H5Dclose(dataset_id);

    // Byn
    ss.str("");ss << "/fields/By/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Byn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Bzn
    ss.str("");ss << "/fields/Bz/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bzn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ex
    ss.str("");ss << "/fields/Ex/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ex[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ey
    ss.str("");ss << "/fields/Ey/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ey[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // Ez
    ss.str("");ss << "/fields/Ez/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ez[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // open the charge density for species
    for (int is = 0; is < ns; is++) {
      ss.str("");ss << "/moments/species_" << is << "/rho/cycle_" << lastcycle;
      dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
      status = H5Dclose(dataset_id);
      array4_double& rhons = *rhons_;
      k = 0;
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int jj = 1; jj < nzn - 1; jj++)
            rhons[is][i][j][jj] = temp_storage[k++];
    }

    // close the hdf file
    status = H5Fclose(file_id);
    delete[]temp_storage;
#endif

}

// extracted from Particles3Dcomm.cpp
//
void Collective::read_particles_restart(
    const VCtopology3D* vct,
    int species_number,
    vector_double& u,
    vector_double& v,
    vector_double& w,
    vector_double& q,
    vector_double& x,
    vector_double& y,
    vector_double& z,
    vector_double& t)const
{
#ifdef NO_HDF5
  eprintf("Require HDF5 to read from restart file.");
#else
    if (vct->getCartesian_rank() == 0 && species_number == 0)
    {
      printf("LOADING PARTICLES FROM RESTART FILE in %s/restart.hdf\n",
        getRestartDirName().c_str());
    }
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = getRestartDirName() + "/restart" + ss.str() + ".hdf";
    // hdf stuff
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[1];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      eprintf("couldn't open file: %s\n"
        "\tRESTART NOT POSSIBLE", name_file.c_str());
      //cout << "couldn't open file: " << name_file << endl;
      //cout << "RESTART NOT POSSIBLE" << endl;
    }


    //find the last cycle
    int lastcycle=0;
    dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &lastcycle);
    status = H5Dclose(dataset_id);

    stringstream species_name;
    species_name << species_number;

    ss.str("");ss << "/particles/species_" << species_number << "/x/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id); /* dataspace handle */
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    // get how many particles there are on this processor for this species
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    const int nop = dims_out[0]; // number of particles in this process
    //Particles3Dcomm::resize_SoA(nop);
    {
      //
      // allocate space for particles including padding
      //
      const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
      u.reserve(padded_nop);
      v.reserve(padded_nop);
      w.reserve(padded_nop);
      q.reserve(padded_nop);
      x.reserve(padded_nop);
      y.reserve(padded_nop);
      z.reserve(padded_nop);
      t.reserve(padded_nop);
      //
      // define size of particle data
      //
      u.resize(nop);
      v.resize(nop);
      w.resize(nop);
      q.resize(nop);
      x.resize(nop);
      y.resize(nop);
      z.resize(nop);
      t.resize(nop);
    }
    // get x
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x[0]);
    // close the data set
    status = H5Dclose(dataset_id);

    // get y
    ss.str("");ss << "/particles/species_" << species_number << "/y/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y[0]);
    status = H5Dclose(dataset_id);

    // get z
    ss.str("");ss << "/particles/species_" << species_number << "/z/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z[0]);
    status = H5Dclose(dataset_id);

    // get u
    ss.str("");ss << "/particles/species_" << species_number << "/u/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    status = H5Dclose(dataset_id);

    // get v
    ss.str("");ss << "/particles/species_" << species_number << "/v/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    status = H5Dclose(dataset_id);

    // get w
    ss.str("");ss << "/particles/species_" << species_number << "/w/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &w[0]);
    status = H5Dclose(dataset_id);

    // get q
    ss.str("");ss << "/particles/species_" << species_number << "/q/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &q[0]);

    //if ID is not saved, read in q as ID
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t[0]);

    status = H5Dclose(dataset_id);

    /* get ID
		ss.str("");ss << "/particles/species_" << species_number << "/ID/cycle_" << lastcycle;
		dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t[0]);
		status = H5Dclose(dataset_id);
    */

    status = H5Fclose(file_id);
#endif
}



/*! constructor */
Collective::Collective(int argc, char **argv) {
  if (argc < 2) {
    inputfile = "inputfile";
    RESTART1 = false;
  }
  else if (argc < 3) {
    inputfile = argv[1];
    RESTART1 = false;
  }
  else {
    if (strcmp(argv[1], "restart") == 0) {
      inputfile = argv[2];
      RESTART1 = true;
    }
    else if (strcmp(argv[2], "restart") == 0) {
      inputfile = argv[1];
      RESTART1 = true;
    }
    else {
      cout << "Error: syntax error in mpirun arguments. Did you mean to 'restart' ?" << endl;
      return;
    }
  }
  ReadInput(inputfile);
  init_derived_parameters();
}

void Collective::init_derived_parameters()
{
  /*! fourpi = 4 greek pi */
  fourpi = 16.0 * atan(1.0);
  /*! dx = space step - X direction */
  dx = Lx / (double) nxc;
  /*! dy = space step - Y direction */
  dy = Ly / (double) nyc;
  /*! dz = space step - Z direction */
  dz = Lz / (double) nzc;
  /*! npcel = number of particles per cell */
  npcel = new int[ns+nstestpart];
  /*! np = number of particles of different species */
  //np = new int[ns];
  /*! npMax = maximum number of particles of different species */
  //npMax = new int[ns];

  /* quantities per process */

  // check that procs divides grid
  // (this restriction should be removed).
  //
  if(0==MPIdata::get_rank())
  {
    fflush(stdout);
    bool xerror = false;
    bool yerror = false;
    bool zerror = false;
    if(nxc % XLEN) xerror=true;
    if(nyc % YLEN) yerror=true;
    if(nzc % ZLEN) zerror=true;
    if(xerror) eprintf("XLEN=%d does not divide nxc=%d\n", XLEN,nxc);
    if(yerror) eprintf("YLEN=%d does not divide nyc=%d\n", YLEN,nyc);
    if(zerror) eprintf("ZLEN=%d does not divide nzc=%d\n", ZLEN,nzc);
    fflush(stdout);
    bool error = (xerror||yerror||zerror);
    // Comment out this check if your postprocessing code does not
    // require the field output subarrays to be the same size.
    // Alternatively, you could modify the output routine to pad
    // with zeros...
    //if(error)
    //{
    //  eprintf("For WriteMethod=default processor dimensions "
    //          "must divide mesh cell dimensions");
    //}
  }

  int num_cells_r = nxc*nyc*nzc;
  //num_procs = XLEN*YLEN*ZLEN;
  //ncells_rs = nxc_rs*nyc_rs*nzc_rs;

  for (int i = 0; i < (ns+nstestpart); i++)
  {
    npcel[i] = npcelx[i] * npcely[i] * npcelz[i];
    //np[i] = npcel[i] * num_cells;
    //nop_rs[i] = npcel[i] * ncells_rs;
    //maxnop_rs[i] = NpMaxNpRatio * nop_rs[i];
    // INT_MAX is about 2 billion, surely enough
    // to index the particles in a single MPI process:
    //assert_le(NpMaxNpRatio * npcel[i] * ncells_proper_per_proc , double(INT_MAX));
    //double npMaxi = (NpMaxNpRatio * np[i]);
    //npMax[i] = (int) npMaxi;
  }
}

/*! destructor */
Collective::~Collective() {
  //delete[]np;
  delete[]npcel;
  delete[]npcelx;
  delete[]npcely;
  delete[]npcelz;
  //delete[]npMax;
  delete[]qom;

  delete[]uth;
  delete[]vth;
  delete[]wth;

  delete[]u0;
  delete[]v0;
  delete[]w0;

  //delete[]TrackParticleID;

  delete[]rhoINIT;

  delete[]pitch_angle;
  delete[]energy;
}
/*! Print Simulation Parameters */
void Collective::Print() {
  cout << endl;
  cout << "Simulation Parameters" << endl;
  cout << "---------------------" << endl;
  cout << "Number of species    = " << ns << endl;
  for (int i = 0; i < ns; i++)
    cout << "qom[" << i << "] = " << qom[i] << endl;
  cout << "x-Length [di]            = " << Lx << endl;
  cout << "y-Length [di]            = " << Ly << endl;
  cout << "z-Length [di]            = " << Lz << endl;
  cout << "Number of cells (x)      = " << nxc << endl;
  cout << "Number of cells (y)      = " << nyc << endl;
  cout << "Number of cells (z)      = " << nzc << endl;
  cout << "Time step                = " << dt << endl;
  cout << "Number of cycles         = " << ncycles << endl;
  cout << "Results saved in  : " << SaveDirName << endl;
  cout << "Case type         : " << Case << endl;
  cout << "Simulation name   : " << SimName << endl;
  cout << "---------------------" << endl;
  cout << "Collision Parameters" << endl;
  cout << "---------------------" << endl;
  cout << "Include Collisions option: " << collisionProcesses << endl;
  cout << "Collision Cross Section  = " << xSec << endl;
  cout << "Species index for secondary electrons: " << iSecElec << endl;
  cout << "Species index for secondary ions:      " << iSecIon << endl;
  if (ns_pl>0)
  {
    cout << "---------------------" << endl;
    cout << "Neutral Gas Parameters" << endl;
    cout << "---------------------" << endl;
    for (int i = 0; i < ns_pl; i++)
    {
      cout << "nSurf[" << i << "] = " << nSurf[i] << ".   ";
      cout << "hExo[" << i << "] = "  << hExo[i] <<  ".   ";
      cout << "fExo[" << i << "] = "  << fExo[i] <<  ".   ";
      cout << endl;
    }
  } 
  cout << "---------------------" << endl;
  cout << "Check Simulation Constraints" << endl;
  cout << "---------------------" << endl;
  cout << "Accuracy Constraint:  " << endl;
  for (int i = 0; i < 1; i++) {
    cout << "u_th < dx/dt species " << i << ".....";
    if (uth[i] < (dx / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;

    cout << "v_th < dy/dt species " << i << "......";
    if (vth[i] < (dy / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;

    cout << "w_th < dz/dt species " << i << "......";
    if (wth[i] < (dz / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;  }
  cout << endl;
  cout << "Finite Grid Stability Constraint:  ";
  cout << endl;
  for (int is = 0; is < 1; is++) {
    if (uth[is] * dt / dx > .1)
      cout << "OK u_th*dt/dx (species " << is << ") = " << uth[is] * dt / dx << " > .1" << endl;
    else
      cout << "WARNING. u_th*dt/dx (species " << is << ") = " << uth[is] * dt / dx << " < .1" << endl;

    if (vth[is] * dt / dy > .1)
      cout << "OK v_th*dt/dy (species " << is << ") = " << vth[is] * dt / dy << " > .1" << endl;
    else
      cout << "WARNING. v_th*dt/dy (species " << is << ") = " << vth[is] * dt / dy << " < .1"  << endl;

    if (wth[is] * dt / dz > .1)
      cout << "OK w_th*dt/dz (species " << is << ") = " << wth[is] * dt / dz << " > .1" << endl;
    else
      cout << "WARNING. w_th*dt/dz (species " << is << ") = " << wth[is] * dt / dz << " < .1"  << endl;
  }

  cout << endl;
  cout << "Collision Probability Constraint:  ";
  cout << endl;
  // Need to account for multiple electron species
  //  where qom < 0
  // Need some parameter for the neutral density to be included to derive an expression for the mfp.
  // For inhomogeneous simulation, where do we want to use for neutral density??
  // Require that vth * dt << mfp = 1/(n * sigma) -> What if n varies by 2 orders of magnitude?
  for (int is = 0; is < ns; is++) 
  {
    if (qom[is] < 0)
    {
      double vthMag2 = uth[is]*uth[is] + vth[is]*vth[is] + wth[is]*wth[is];
      double vthMag = sqrt(vthMag2);
      double totMfpInv = 0.0; // Inverse of the mean free path i.e Sum( n * xSec)

      // Sum over collisional neutral species
      for (int iNeut = 0; iNeut < ns_pl; iNeut ++) 
        totMfpInv += nSurf[iNeut] * xSec;
  
      if ((vthMag * dt * totMfpInv < .1)) // If collisions are unlikely to occur
        cout << "OK |u_th| * dt / mfp (species " << is << ") = " << vthMag * dt * totMfpInv << " < .1" << endl;
      else
        cout << "WARNING. |u_th| * dt / mfp (species " << is << ") = " << vthMag * dt * totMfpInv << " > .1" << endl;
    }
  }
  

  cout << "\n" << "n_layers_sal: " << n_layers_sal << "\n";
  cout << "\n" << "Nx/XLEN: " << nxc/XLEN << "\n";
  cout << "\n" << "Ny/YLEN: " << nyc/YLEN << "\n";
  cout << "\n" << "Nz/ZLEN: " << nzc/ZLEN << "\n";

  if (yes_sal){
    cout << "\nREQUIRE SAL_length/(Vthi*dt)>>1 = " << n_layers_sal*dz/dt/vth[1] << endl;
    if (n_layers_sal>=(nxc/XLEN)) 
      eprintf("ERROR: n_layers_sal bigger than Nx/XLEN")
    if (n_layers_sal>=(nyc/YLEN)) 
      eprintf("ERROR: n_layers_sal bigger than Ny/YLEN")
     if (n_layers_sal>=(nzc/ZLEN)) 
      eprintf("ERROR: n_layers_sal bigger than Nz/ZLEN")
   }


}
/*! Print Simulation Parameters in tree.xml */
// SPASE-compliant metadata in xml file
void Collective::save() {

  // path to outputfile
  string temp;  
  temp = SaveDirName + "/tree.xml";
  ofstream my_file(temp.c_str());

  // day of the simulation
  time_t t = time(0); 
  tm* now = localtime(&t);
  stringstream year, month, day, hour, min, sec;
  year << setw(4) << setfill('0') << (now->tm_year+1900);
  month << setw(2) << setfill('0') << (now->tm_mon+1);
  day << setw(2) << setfill('0') << now->tm_mday;
  hour << setw(2) << setfill('0') << now->tm_hour;
  min << setw(2) << setfill('0') << now->tm_min;
  sec << setw(2) << setfill('0') << now->tm_sec;

  // string for later
  string name_sp, name_ext;

  // url to XML standard
  string url_standard = "http://www.w3.org/2001/XMLSchema-instance";
  // url to spase website - used to check scheme tree
  string url_spase = "http://www.spase-group.org/data/schema";
  // url to server where data are stored - TO BE UPDATED
  string url_server = "http://impex.latmos.ipsl.fr";
  // date of release of the code
  string release_date = "2024-01-01T00:00:00.000";
  // version of the code release
  string vers = "0.0.1";
  // date of start of this simulation
  string date_simu = year.str()+"-"+month.str()+"-"+day.str()+"T"+hour.str()+":"+min.str()+":"+sec.str()+".000";
  // name of this simulation
  string name = "Mercury_"+day.str()+"_"+month.str()+"_"+year.str();
  // name of the ouput file (for now only mag field at t=0)
  string name_output = SimName + "_" + "B" + "_" + "0"; 

  // Header of the xml file - should NOT be changed
  my_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  my_file << "<Spase>" << endl; 
  my_file << "\t<Version>"<< vers <<"</Version>" << endl;

  // Summary information on the simu Model
  my_file << "\t<SimulationModel>" << endl;
    my_file << "\t\t<ResourceID>spase://IMPEX/SimulationModel/Lagrange</ResourceID>" << endl;
    my_file << "\t\t<ResourceHeader>" << endl;
      my_file << "\t\t\t<ResourceName>PIC_Lagrange</ResourceName>" << endl;
      my_file << "\t\t\t<ReleaseDate>" << release_date << "</ReleaseDate>" << endl;
      my_file << "\t\t\t<Description>" << endl;
        my_file << "\t\t\t\tPIC simulation model developped at OCA for plasma interaction with celestial neutral environments (planets or moons)" << endl;
      my_file << "\t\t\t</Description>" << endl;
      my_file << "\t\t\t<Contact>" << endl;
        my_file << "\t\t\t\t<PersonID>Lagrange</PersonID>" << endl;
        my_file << "\t\t\t\t<Role>DataProducer</Role>" << endl;
      my_file << "\t\t\t</Contact>" << endl;
      my_file << "\t\t\t<InformationURL>" << endl;
        my_file << "\t\t\t\t<URL>" << url_server << "</URL>" << endl;
      my_file << "\t\t\t</InformationURL>" << endl;
    my_file << "\t\t</ResourceHeader>" << endl;
    my_file << "\t\t<SimulationType>PIC</SimulationType>" << endl;
    my_file << "\t\t<CodeLanguage>C/C++</CodeLanguage>" << endl;
  my_file << "\t</SimulationModel>" << endl;

  // Information of the rep where the simulations are stored
  my_file << "\t<Repository>" << endl; 
    my_file << "\t\t<ResourceID>spase://IMPEX/Repository/Lagrange</ResourceID>" << endl; 
    my_file << "\t\t<ResourceHeader>" << endl; 
      my_file << "\t\t\t<ResourceName>Lagrange_PIC_Simulation_Database</ResourceName>" << endl; 
      my_file << "\t\t\t<ReleaseDate>" << release_date << "</ReleaseDate>" << endl;
      my_file << "\t\t\t<Description/>" << endl;
      my_file << "\t\t\t<Contact>" << endl; 
        my_file << "\t\t\t\t<PersonID>Lagrange</PersonID>" << endl; 
        my_file << "\t\t\t\t<Role>DataProducer</Role>" << endl; 
      my_file << "\t\t\t</Contact>" << endl; 
    my_file << "\t\t</ResourceHeader>" << endl;
    my_file << "\t\t<AccessURL>" << endl; 
      my_file << "\t\t\t<URL>" << url_server << "</URL>" << endl; 
    my_file << "\t\t</AccessURL>" << endl; 
  my_file << "\t</Repository>" << endl; 

  // Information on the run we are doing now
  my_file << "\t<SimulationRun>" << endl;
    my_file << "\t\t<ResourceID>spase://IMPEX/SimulationRun/Lagrange/" << name << "</ResourceID>" << endl;
    my_file << "\t\t<ResourceHeader>" << endl;
      my_file << "\t\t\t<ResourceName>" << name << "</ResourceName>" << endl;
      my_file << "\t\t\t<ReleaseDate>" << date_simu << "</ReleaseDate>" << endl;
      my_file << "\t\t\t<Description/>" << endl;
      my_file << "\t\t\t<Contact>" << endl;
        my_file << "\t\t\t\t<PersonID>Lagrange</PersonID>" << endl;
        my_file << "\t\t\t\t<Role>DataProducer</Role>" << endl;
      my_file << "\t\t\t</Contact>" << endl;
    my_file << "\t\t</ResourceHeader>" << endl;
    my_file << "\t\t<Model>" << endl;
      my_file << "\t\t\t<ModelID>spase://IMPEX/SimulationModel/Lagrange</ModelID>" << endl;
    my_file << "\t\t</Model>" << endl;

    // info on time of simulation --> TO BE CHANGED with real values
    my_file << "\t\t<TemporalDependence>Yes</TemporalDependence>" << endl;
    my_file << "\t\t<SimulatedRegion>Mercury</SimulatedRegion>" << endl;
    my_file << "\t\t<LikelihoodRating>Strong</LikelihoodRating>" << endl;
    my_file << "\t\t<SimulationTime>" << endl;
      my_file << "\t\t\t<Duration>PT2072.854S</Duration>" << endl;  //put real tim of the simu (numerical)? or simu world time?
      my_file << "\t\t\t<TimeStart>00:00:00</TimeStart>" << endl;
      my_file << "\t\t\t<TimeStop>00:00:00</TimeStop>" << endl;
      my_file << "\t\t\t<TimeStep>PT0.115S</TimeStep>" << endl;
    my_file << "\t\t</SimulationTime>" << endl;

    // info of simu box
    my_file << "\t\t<SimulationDomain>" << endl;
      my_file << "\t\t\t<CoordinateSystem>" << endl;
        my_file << "\t\t\t\t<CoordinateRepresentation>Cartesian</CoordinateRepresentation>" << endl;
        my_file << "\t\t\t\t<CoordinateSystemName>Box</CoordinateSystemName>" << endl;
      my_file << "\t\t\t</CoordinateSystem>" << endl;
      my_file << "\t\t\t<SpatialDimension>3</SpatialDimension>" << endl;
      my_file << "\t\t\t<FieldDimension>3</FieldDimension>" << endl;
      my_file << "\t\t\t<Units>di</Units>" << endl;
      my_file << "\t\t\t<UnitsConversion> " << 1./L_square << " &gt; R </UnitsConversion>" << endl;
      my_file << "\t\t\t<CoordinatesLabel>X Y Z</CoordinatesLabel>" << endl;
      my_file << "\t\t\t<ValidMin>0  0  0</ValidMin>" << endl;
      my_file << "\t\t\t<ValidMax>" << Lx << " " << Ly << " " << Lz << "</ValidMax>" << endl;
      my_file << "\t\t\t<GridStructure>Constant</GridStructure>" << endl;
      my_file << "\t\t\t<GridCellSize>" << Lx/nxc << " " << Ly/nyc << " " << Lz/nzc << "</GridCellSize>" << endl;
      my_file << "\t\t\t<Symmetry>Plane</Symmetry>" << endl;
      my_file << "\t\t\t<BoundaryConditions>" << endl;
        my_file << "\t\t\t\t<ParticleBoundary>" << endl;
          my_file << "\t\t\t\t\t<FrontWall> absorbing </FrontWall>" << endl; // definisci meglio qui le BC
          my_file << "\t\t\t\t\t<BackWall> absorbing </BackWall>" << endl;
          my_file << "\t\t\t\t\t<SideWall> absorbing </SideWall>" << endl;
          my_file << "\t\t\t\t\t<Obstacle> absorbing </Obstacle>" << endl;
        my_file << "\t\t\t\t</ParticleBoundary>" << endl;
        my_file << "\t\t\t\t<FieldBoundary>" << endl;
          my_file << "\t\t\t\t\t<FrontWall> IMF </FrontWall>" << endl;		// anche queste... Le vogliamo fissare o lasciamo qualche scelta allo user?
          my_file << "\t\t\t\t\t<BackWall> Neuman zero-gradient </BackWall>" << endl;
          my_file << "\t\t\t\t\t<SideWall> periodic </SideWall>" << endl;
          my_file << "\t\t\t\t\t<Obstacle> absorbing </Obstacle>" << endl;
        my_file << "\t\t\t\t</FieldBoundary>" << endl;
      my_file << "\t\t\t</BoundaryConditions>" << endl;
    my_file << "\t\t</SimulationDomain>" << endl;

    // Information on the planet / body simulated
    // 2- altitude exobas add real value and nsurf
    // 4- add dipole offset
    // 5- add info on dipole field (B1, harmonic expansion?)
    my_file << "\t\t<RegionParameter>" << endl;
      my_file << "\t\t\t<SimulatedRegion>Mercury</SimulatedRegion>" << endl;
      my_file << "\t\t\t<Radius Units=di>" << L_square << "</Radius>" << endl;
      my_file << "\t\t\t<Property>" << endl;
        my_file << "\t\t\t\t<Name>Planet center</Name>" << endl;
        my_file << "\t\t\t\t<Description>Coordinates of planet center w.r.t. box reference frame</Description>" << endl;
        my_file << "\t\t\t\t<PropertyQuantity>Cartesion Coordinates</PropertyQuantity>" << endl;
        my_file << "\t\t\t<Units>di</Units>" << endl;
        my_file << "\t\t\t\t<PropertyValue>" << x_center << " " << y_center << " " << z_center << "</PropertyValue>" << endl;
      my_file << "\t\t\t</Property>" << endl;
    my_file << "\t\t</RegionParameter>" << endl;

    // Information on input field B_SW
    my_file << "\t\t<InputField>" << endl;
      my_file << "\t\t\t<Name>Magnetic Field SW</Name>" << endl;
      my_file << "\t\t\t<Description>Interplanetary Magnetic Field</Description>" << endl;
      my_file << "\t\t\t<SimulatedRegion>Heliosphere</SimulatedRegion>" << endl;
      my_file << "\t\t\t<FieldQuantity>Magnetic</FieldQuantity>" << endl;
      my_file << "\t\t\t<Units>1/wpi</Units>" << endl;
      my_file << "\t\t\t<InputLabel>Bx By Bz</InputLabel>" << endl;
      my_file << "\t\t\t<FieldValue>" << B0x << " " << B0y << " " << B0z << "</FieldValue>" << endl;
    my_file << "\t\t</InputField>" << endl;

    // Information on input field E_SW
    my_file << "\t\t<InputField>" << endl;
      my_file << "\t\t\t<Name>Electric Field SW</Name>" << endl;
      my_file << "\t\t\t<Description>Interplanetary Motional Electric Field</Description>" << endl;
      my_file << "\t\t\t<SimulatedRegion>Heliosphere</SimulatedRegion>" << endl;
      my_file << "\t\t\t<FieldQuantity>Electric</FieldQuantity>" << endl;
      my_file << "\t\t\t<Units>1/wpi</Units>" << endl;
      my_file << "\t\t\t<InputLabel>Ex Ey Ez</InputLabel>" << endl;
      my_file << "\t\t\t<FieldValue>" << E0x << " " << E0y << " " << E0z << "</FieldValue>" << endl;
    my_file << "\t\t</InputField>" << endl;

    // Information on input field B_planet
    my_file << "\t\t<InputField>" << endl;
      my_file << "\t\t\t<Name>Magnetic Field Planet</Name>" << endl;
      my_file << "\t\t\t<Description>Planet Intrinsic Magnetic Field</Description>" << endl;
      my_file << "\t\t\t<SimulatedRegion>Magnetosphere</SimulatedRegion>" << endl;
      my_file << "\t\t\t<FieldQuantity>Magnetic</FieldQuantity>" << endl;
      my_file << "\t\t\t<Units>1/wpi</Units>" << endl;
      my_file << "\t\t\t<InputLabel>Bz</InputLabel>" << endl;
      my_file << "\t\t\t<FieldValue>" << B1z << "</FieldValue>" << endl;
    my_file << "\t\t</InputField>" << endl;

    // Information on input Particles (both SW and planet)
    for (int is=0; is<ns; is++)
    {
      if (qom[is]>0)
      {
        name_sp = "i";
        name_ext = "ions";
      }
      if (qom[is]<0)
      {
        name_sp = "e";
        name_ext = "electrons";
      }
      my_file << "\t\t<InputPopulation>" << endl;
        my_file << "\t\t\t<Name>Species " << name_sp << is << "</Name>" << endl;
        my_file << "\t\t\t<ParticleType>" << name_ext << "</ParticleType>" << endl;
        my_file << "\t\t\t<PopulationMassNumber>" << fabs(1./qom[is]) <<"</PopulationMassNumber>" << endl;
        my_file << "\t\t\t<PopulationChargeState>" << (-2*signbit(qom[is]))+1 << "</PopulationChargeState>" << endl;   // assume charge state is always +-1
        my_file << "\t\t\t<PopulationTemperature Units=mic2>" << fabs(1./qom[is])*(uth[is]*uth[is]+vth[is]*vth[is]+wth[is]*wth[is])<< "</PopulationTemperature>" << endl; // assume No anisotropy in the solar wind pcls (can be changed...)
        my_file << "\t\t\t<PopulationFlowSpeed Units=c>" << u0[is] << " " << v0[is] << " " << w0[is] << "</PopulationFlowSpeed>" << endl;
        my_file << "\t\t\t<Distribution>Maxwellian</Distribution>" << endl;
        if (is<ns_sw)
        my_file << "\t\t\t	<PopulationDensity Units=nsw>" << rhoINIT[is] <<  "</PopulationDensity>" << endl;
	else
	{
        my_file << "\t\t\t<Property>" << endl;
          my_file << "\t\t\t\t<Name>Surface Density</Name>" << endl;
          my_file << "\t\t\t\t<Description>Surface Density of neutrals</Description>" << endl;
          my_file << "\t\t\t\t<PropertyQuantity>Density</PropertyQuantity>" << endl;
          my_file << "\t\t\t\t<Units>nsw</Units>" << endl;
          my_file << "\t\t\t\t<PropertyValue>" << nSurf[is-ns_sw] << "</PropertyValue>" << endl;
        my_file << "\t\t\t</Property>" << endl;
        my_file << "\t\t\t<Property>" << endl;
          my_file << "\t\t\t\t<Name>Exobase Altitude</Name>" << endl;
          my_file << "\t\t\t\t<Description>Altitude of the exobase</Description>" << endl;
          my_file << "\t\t\t\t<PropertyQuantity>Positional</PropertyQuantity>" << endl;
          my_file << "\t\t\t\t<Units>di</Units>" << endl;
          my_file << "\t\t\t\t<PropertyValue>" << hExo[is-ns_sw] << "</PropertyValue>" << endl;
        my_file << "\t\t\t</Property>" << endl;
	}	
      }
      my_file << "\t\t</InputPopulation>" << endl;

    // Info on electron-neutral collision
    for (int icoll=nCollProcesses-1; icoll>=0; icoll-=1)
    {
    if (icoll==nCollProcesses-1)
      Emax = 1.0;
    else
      Emax = E_th_el[icoll+1];
    Emin = E_th_el[icoll];
    my_file << "\t\t<InputProcess>" << endl;
      my_file << "\t\t\t<Name>Electron-Neutral Collision #" << icoll << "</Name>" << endl;
      my_file << "\t\t\t<Description> Impact of electrons with energy in range " << Emin << " " << Emax << " [mic2]</Description>" << endl;
      my_file << "\t\t\t<ProcessType>Collision</ProcessType>" << endl;
      my_file << "\t\t\t<ProcessModel>Pete2022</ProcessModel>" << endl;  // look for a reference for this collision module (if existing ???)
      my_file << "\t\t\t<Units>1/nsw/di</Units>" << endl;
      my_file << "\t\t\t<ProcessCoefficient>" << xSec << "</ProcessCoefficient>" << endl;
      my_file << "\t\t\t<ProcessCoeffType>Cross Section</ProcessCoeffType>" << endl;
    my_file << "\t\t</InputProcess>" << endl;
    }

    // Info on plasma creation via photoionization
    for (int iph=ns_sw; iph<(ns_sw+ns_pl); iph++)
    {
    my_file << "\t\t<InputProcess>" << endl;
      my_file << "\t\t\t<Name>Creation of plasma of species " << iph << " via photoionization</Name>" << endl;
      my_file << "\t\t\t<ProcessType>Photoionization</ProcessType>" << endl;
      my_file << "\t\t\t<ProcessModel>Plasma creation with constant rate.</ProcessModel>" << endl;
      my_file << "\t\t\t<Units>1/wpi</Units>" << endl;
      my_file << "\t\t\t<ProcessCoefficient>" << fExo[iph-ns_sw] << "</ProcessCoefficient>" << endl;
      my_file << "\t\t\t<ProcessCoeffType>frequency</ProcessCoeffType>" << endl;
    my_file << "\t\t</InputProcess>" << endl;
    }

    // Derived parameters : Alfven speed, plasma beta, mach number
    my_file << "\t\t<InputParameter>" << endl;
      my_file << "\t\t\t<Name>Derived Parameters</Name>" << endl;
      my_file << "\t\t\t<ParameterQuantity>Other</ParameterQuantity>" << endl;
      my_file << "\t\t\t<Property>" << endl;
        my_file << "\t\t\t\t<Name>Alfven Speed</Name>" << endl;
        my_file << "\t\t\t\t<PropertyQuantity>AlfvenVelocity</PropertyQuantity>" << endl;
        my_file << "\t\t\t\t<Units>c</Units>" << endl;
        my_file << "\t\t\t\t<PropertyValue>" << sqrt(B0x*B0x+B0y*B0y+B0z*B0z) << "</PropertyValue>" << endl;
      my_file << "\t\t\t</Property>" << endl;
      my_file << "\t\t\t<Property>" << endl;
        my_file << "\t\t\t\t<Name>Plasma Beta Ions i1</Name>" << endl;
        my_file << "\t\t\t\t<PropertyQuantity>RatioPlasma-to-MagneticPressureIons</PropertyQuantity>" << endl;
        my_file << "\t\t\t\t<PropertyValue>" << 2.*(uth[1]*uth[1]+vth[1]*vth[1]*+wth[1]*wth[1])/(B0x*B0x+B0y*B0y+B0z*B0z) << "</PropertyValue>" << endl;
      my_file << "\t\t\t</Property>" << endl;
      my_file << "\t\t\t<Property>" << endl;
        my_file << "\t\t\t\t<Name>Alfven Mach Number Ions i1</Name>" << endl;
        my_file << "\t\t\t\t<PropertyQuantity>MachNumber</PropertyQuantity>" << endl;
        my_file << "\t\t\t\t<PropertyValue>" << sqrt((u0[1]*u0[1]+v0[1]*v0[1]+w0[1]*w0[1])/(B0x*B0x+B0y*B0y+B0z*B0z)) << "</PropertyValue>" << endl;
      my_file << "\t\t\t</Property>" << endl;
      my_file << "\t\t\t<Property>" << endl;
        my_file << "\t\t\t\t<Name>Sonic Mach Number Ions i1</Name>" << endl;
        my_file << "\t\t\t\t<PropertyQuantity>MachNumber</PropertyQuantity>" << endl;
        my_file << "\t\t\t\t<PropertyValue>" <<  sqrt((u0[1]*u0[1]+v0[1]*v0[1]+w0[1]*w0[1])/(uth[1]*uth[1]+vth[1]*vth[1]*+wth[1]*wth[1])) << "</PropertyValue>" << endl;
      my_file << "\t\t\t</Property>" << endl;
    my_file << "\t\t</InputParameter>" << endl;

  my_file << "\t</SimulationRun>" << endl; 

  // Example for output 3D magnetic field 
  // need a loop for all quantities printed and shared
  my_file << "\t<NumericalOutput>" << endl;
    my_file << "\t\t<ResourceID>spase://IMPEX/NumericalOutput/Lagrange/" << name << "/Mag/3D</ResourceID>" << endl;
    my_file << "\t\t<ResourceHeader>" << endl;
      my_file << "\t\t\t<ResourceName>Mag/3D</ResourceName>" << endl;
      my_file << "\t\t\t<ReleaseDate>" << date_simu << "</ReleaseDate>" << endl;
      my_file << "\t\t\t<Description/>" << endl;
      my_file << "\t\t\t<Contact>" << endl;
        my_file << "\t\t\t\t<PersonID>Lagrange</PersonID>" << endl;
        my_file << "\t\t\t\t<Role>DataProducer</Role>" << endl;
      my_file << "\t\t\t</Contact>" << endl;
    my_file << "\t\t</ResourceHeader>" << endl;
    my_file << "\t\t<AccessInformation>" << endl;
      my_file << "\t\t\t<RepositoryID>spase://IMPEX/Repository/Lagrange</RepositoryID>" << endl;
      my_file << "\t\t\t<AccessURL>" << endl;
        my_file << "\t\t\t\t<URL>" << url_server << "</URL>" << endl;
      my_file << "\t\t\t</AccessURL>" << endl;
      my_file << "\t\t\t<Format>HDF5</Format>" << endl;
    my_file << "\t\t</AccessInformation>" << endl;
    my_file << "\t\t<MeasurementType>MagneticField</MeasurementType>" << endl;
    my_file << "\t\t<SpatialDescription>" << endl;
      my_file << "\t\t\t<Dimension>3</Dimension>" << endl;
      my_file << "\t\t\t<CoordinateSystem>" << endl;
        my_file << "\t\t\t\t<CoordinateRepresentation>Cartesian</CoordinateRepresentation>" << endl;
        my_file << "\t\t\t\t<CoordinateSystemName>box</CoordinateSystemName>" << endl;
      my_file << "\t\t\t</CoordinateSystem>" << endl;
      my_file << "\t\t\t<Units>di</Units>" << endl;
      my_file << "\t\t\t<UnitsConversion>" << 1./L_square << " &gt; R</UnitsConversion>" << endl;
      my_file << "\t\t\t<RegionBegin> 0. 0. 0. </RegionBegin>" << endl;
      my_file << "\t\t\t<RegionEnd> " << Lx << " " << Ly << " " << Lz << " </RegionEnd>" << endl;  // for now keep box units - to be check if need to add here MSO coord.
    my_file << "\t\t</SpatialDescription>" << endl;
    my_file << "\t\t<SimulatedRegion>Mercury</SimulatedRegion>" << endl;
    my_file << "\t\t<InputResourceID>spase://IMPEX/SimulationRun/Lagrange/" << name << "</InputResourceID>" << endl;
    my_file << "\t\t<Parameter>" << endl;
      my_file << "\t\t\t<Name>MagneticField</Name>" << endl;
      my_file << "\t\t\t<ParameterKey>Bx By Bz</ParameterKey>" << endl;
      my_file << "\t\t\t<Units>1/wpi</Units>" << endl;
      my_file << "\t\t\t<UnitsConversion> " << 1./sqrt(B0x*B0x+B0y*B0y+B0z*B0z) << " &gt; Bsw </UnitsConversion>" << endl; // TBC if need to add change to nT (is this possible?)
      my_file << "\t\t\t<Field>" << endl;
        my_file << "\t\t\t\t<Qualifier>Vector</Qualifier>" << endl;
        my_file << "\t\t\t\t<FieldQuantity>Magnetic</FieldQuantity>" << endl;
      my_file << "\t\t\t</Field>" << endl;
    my_file << "\t\t</Parameter>" << endl;
    my_file << "\t\t<SimulationProduct>3DCubes</SimulationProduct>" << endl;
  my_file << "\t</NumericalOutput>" << endl;

  // create a granule corrsponding to one file printed
  my_file << "\t<Granule>" << endl;
    my_file << "\t\t<ResourceID>spase://IMPEX/Granule/Lagrange/" << name << "/Mag/3D/" << name_output << "</ResourceID>" << endl;
    my_file << "\t\t<ReleaseDate>" << date_simu << "</ReleaseDate>" << endl;
    my_file << "\t\t<ParentID>spase://IMPEX/NumericalOutput/Lagrange/" << name << "/Mag/3D</ParentID>" << endl;
    my_file << "\t\t\t<RegionBegin> 0. 0. 0. </RegionBegin>" << endl;
    my_file << "\t\t\t<RegionEnd> " << Lx << " " << Ly << " " << Lz << " </RegionEnd>" << endl;  // for now keep box units - to be check if need to add here MSO coord.
    my_file << "\t\t<Source>" << endl;
      my_file << "\t\t\t<SourceType>Data</SourceType>" << endl;
      my_file << "\t\t\t<URL>" << url_server << "/" << name << "/" << name_output << ".h5</URL>" << endl;
    my_file << "\t\t</Source>" << endl;
  my_file << "\t</Granule>" << endl;

my_file << "</Spase>" << endl;

 my_file.close();
}

