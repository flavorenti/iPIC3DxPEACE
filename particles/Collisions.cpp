/*******************************************************************************************
  Collisions.h -  Collisional processes for particles.
  -------------------
developers: Peter Stephenson
 ********************************************************************************************/


#include <iostream>
// #include <fstream>
// #include <sstream>
#include <math.h>
#include "Particles3D.h"
#include "Collective.h"
#include "Collisions.h"
#include "ipicmath.h"


Collisions::~Collisions()
{

}

Collisions::Collisions(CollectiveIO * col, VirtualTopology3D *vct, Grid * grid)
{
// Load simulation parameters
Parameters::init_parameters();
// Load cross sections for simulation
xSec = col->getxSec();
iSecElec = col->getiSecElec();
iSecIon = col->getiSecIon();
nCollProcesses = col->getnCollProcesses();
nIoniColls = col->getnIoniColls();
collStepSkip = col->getcollStepSkip();
E_th_el = new double[nCollProcesses];
for (int i = 0; i < nCollProcesses; i++){
  E_th_el[i] = col->getE_th_el(i);
}


dt = col->getDt();
ns = col->getNs();   
// std::cout << "\n Wowwwwww \n \n \n Such Collide \n";

// Set up ionized Pls object
// Allocation of particles
  eImpIoniPls = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int i = 0; i < ns; i++)
  {
    new(&eImpIoniPls[i]) Particles3D(i,col,vct,grid);
  }

} 

// Apply Collisions from different species. 
void Collisions::Collide(int species, Particles3D *part, CollectiveIO * col)
{
  // cout << endl << "Collide for Species: " << species << endl;
  // Get charge of the particle species
  double qom = col->getQOM(species); 
  double vScale = 1; // velocity scaling for ions
  if (qom<0)// If electrons
  {
  int nPls = part[species].getNOP();

  // Scaling of velocity required when using heavier electrons
  // Not needed for ions
  vScale = sqrt(qom_eReal / qom); // Electrons

  // cout << "Number of electrons in Species " << species << ": " << nPls << endl;
  for (int pidx = 0; pidx < nPls; pidx++) {
    // COmpact version
    double vMag = velMagnitude(species, part, pidx); // Find particle velocity
    double pColl = ProbColl(vMag, vScale); // Calculate collision probability
    double rColl =  sample_clopen_u_double();// Random number for coll
    // cout << "pColl: " << pColl << endl;


    // Epl = 0.5 * pow(vMag, 2) / abs(qom); // Change mass to abs(1/QOM)
    // cout << "Vmag: " << vMag << "Electon Energy: " << Epl << endl;


    if (rColl < pColl) // Particle undergoes collision
    {
      // cout << "Pl " << pidx  << ". pColl, rColl: " << pColl << " , " << rColl <<  endl;
      // cout << "Electron Collides \n";
      CollideElectron(species, part, pidx, vMag, col); // Modify electron velocity
    }

    // double upl = part[species].getU(pidx);
    // double vpl = part[species].getV(pidx);
    // double wpl = part[species].getW(pidx);
    // double vMag = sqrt(pow(upl,2) + pow(vpl,2) + pow(wpl,2));
    
    // // Calculate neutral density for each particle
    // // Find particle positions
    // const double xpl = part[species].getX(pidx);
    // const double ypl = part[species].getY(pidx);
    // const double zpl = part[species].getZ(pidx);
    // nNeutral = 1;

    
    // const double tau = nNeutral * xSec * vMag * dt;
    // const double pColl = 1 - exp(- tau);
    
    
    // // std::cout << "upl = " << upl << vpl << wpl <<endl;
    // // std::cout << "xpl = " << xpl << ypl<< zpl << endl;
    // // std::cout << "xSec = " << xSec << endl;
    // // std::cout << "dt = " << nNeutral << endl;
    // // std::cout << "nNeutral = " << nNeutral << endl;
    // // std::cout << "xSec = " << xSec << endl;
    // // std::cout << "vMag = " << vMag << endl;
    // // std::cout << "tau = " << tau << endl;
    // // std::cout << "pColl = " << pColl << endl;

    // // Get random numbers for each particle
    // double rColl =  sample_clopen_u_double();
    // // Implement collisions
    // if (rColl < pColl)
    // {
    //   // CollideElectron(species, part, pidx, vMag);
    //   // cout << "pColl = " << pColl << "\n";
    //   // cout << "rColl = " << rColl << "\n";

    // // // Find particle Energy and energy loss for each colliding particle
    // Epl = 0.5 * mPl * pow(vMag, 2);
    // if (Epl > 13.5) Eth = 13.5 + 0.5 * (Epl - 13.5);//Ionization
    // else if (Epl > 8) Eth = 8; // Electronic Transitions
    // else if (Epl > 0.4) Eth = 0.4; //Vibrational Collisions
    // else Eth = 0; 
    // // // Reduce particle energies 
    // Epl_new = Epl - Eth;
    // // // Determine new particle velocities
    // double vFact = sqrt(Epl_new /Epl);
    // // cout << "Old Energy:" << Epl << "\n";
    // // cout << "New Energy:" << Epl_new << "\n";
    // // cout << "Old Velocities:" << upl << " " << vpl << " " << wpl << "\n ";
    // upl *= vFact;
    // vpl *= vFact;
    // wpl *= vFact;
    // // cout << "New Velocities:" << upl << " " << vpl << " " << wpl << "\n ";
    // // Set particle velocities
    // part[species].setU(pidx, upl);
    // part[species].setV(pidx, vpl);
    // part[species].setW(pidx, wpl);
    // // cout << "Particle underwent collision!!! \n";
    // }
  }
  }
}


// Modififes electron velocities in case of collision. 
// Also saves properties of ionized pls. 
void Collisions::CollideElectron(int species, Particles3D *part, int pidx, double vMag, CollectiveIO * col)
{
 // Find particle Energy and energy loss for each colliding particle
 // Energy values need to be rescaled into code units. 
  double qom = col->getQOM(species);
  Epl = 0.5 * pow(vMag, 2) / abs(qom); // Change mass to abs(1/QOM)

  // E_th_el>0
  // Run through collisions with decreasing threshold energy
  isIoni = false;
  for (int iColl = 0; iColl < nCollProcesses; iColl++)
  {
    if (Epl > E_th_el[iColl])
    {
      
      Eth = E_th_el[iColl];
      
      if (iColl < nIoniColls)
      {
        // Splitting energy 50-50 between electrons
        Eth = Eth + 0.5 * (Epl - Eth);
        isIoni = true; 
      }
      cout << "iColl: " << iColl << ". Eth = " << Eth << ". isIoni = " << isIoni << endl;
      break;
      
    }
  }
  // }

//   // cout << "Vmag: " << vMag << "Electon Energy: " << Epl << endl; 
//   if (Epl > 13.5) Eth = 13.5 + 0.5 * (Epl - 13.5);//Ionization - using 50-50 split
//   else if (Epl > 8) Eth = 8; // Electronic Transitions
//   else if (Epl > 0.4) Eth = 0.4; //Vibrational Collisions
//   else Eth = 0; // No energy loss below 0.4 eV
//  // Could add continuous cooling processes for v low energies.

  // Reduce particle energies 
  Epl_new = Epl - Eth;
  // Determine new particle velocities
  double vFact = sqrt(Epl_new /Epl);
  // cout << "OldEnergy:" << Epl << " ";
  // cout << "New Energy:" << Epl_new << "\n";
  // cout << "vFact: " << vFact << "\n";
  // cout << "Old Velocities:" << upl << " " << vpl << " " << wpl << "\n ";
  // cout << "Old Vel:" << upl << " " << vpl << " " << wpl << "\n";
  
  // Scales electron velocities after energy loss
  upl *= vFact;
  vpl *= vFact;
  wpl *= vFact;
  // cout << "New Vel:" << upl << " " << vpl << " " << wpl << "\n";
  // cout << "New Velocities:" << upl << " " << vpl << " " << wpl << "\n ";
  // Set particle velocities
  part[species].setU(pidx, upl);
  part[species].setV(pidx, vpl);
  part[species].setW(pidx, wpl);

// If Ionization collisions
  if (isIoni) 
  {
    // Save properties of new ionized pls.
    // Information from primary electrons
    recordIonizedParticles(species, part, pidx, col);
  }
}

/* Calculate Collision probability in a given step */
double Collisions::ProbColl(double vMag, double vScale)
{
  // nNeutral - neutral density of water (?) 
  // Should be same as used for photoionization module

  // xSec - Cross section of collision. Set in input file
  // vScale - Scaling ratio to correct for heavy electrons
  double tau = nNeutral * xSec * vScale * vMag * dt;

  /* If steps skipped then need to scale up collision probability
  by an appropriate amount */
  tau = tau * collStepSkip;
  double pColl = 1 - exp(-tau);
  return pColl;
}

/* Returns the velocity magnitude for a particle */
double Collisions::velMagnitude(int species, Particles3D *part, int pidx)
{
  // vMag = velMagnitude(part[species], pidx);
        
  // Find particle velocity
  upl = part[species].getU(pidx);
  vpl = part[species].getV(pidx);
  wpl = part[species].getW(pidx);
  double vMag = sqrt(pow(upl,2) + pow(vpl,2) + pow(wpl,2));

  return vMag;
}

void Collisions::ImpactIonization()
{
  // Save Pl properties after ionization

}

/* Records the propertpies of secondary electrons and ions in case of ionization */
void Collisions::recordIonizedParticles(int species, Particles3D *part, int pidx, CollectiveIO * col)
{
  // Save Pl properties after ionization
  /* Need to retain - 
  position and energy/velocity for secondary electrons 
  Position only for ions - created stationary?
  */
  // species - index of the ionizing species.
  /* Get information about the ionizing particle*/
  // cout << "\n Saving Ionized Pls Info \n";

  // Information input from primary electrons.

  // Position and velocity of ionizing particle.
  double x = part[species].getX(pidx);
  double y = part[species].getY(pidx);
  double z = part[species].getZ(pidx);
  
  double u = part[species].getU(pidx);
  double v = part[species].getV(pidx);
  double w = part[species].getW(pidx);
  // cout << "x: " << x << " " << y << " " << z << "\n";
  // cout << "v: " << u << " " << v << " " << w << "\n";
  // cout << "Pl info retrieved \n";
  

  // Need to remove hard coding here. Retrieve from input file
  int iComEl = 0;
  int iComIon = 1;
  /* Save ionospheric eletrons produced by collisions*/
  // cout << "\n Retrieving charge from collective \n";
  // cout << "iComEl: " << iComEl <<"\n"; 
  // int ns = col->getNs();
  // int nsTestPart = col->getNsTestPart();
  // cout << "Col ns: " << ns << "\n";
  // cout << "Col nsTestPart: " << nsTestPart << "\n";
  

  // record secondary electron
  // double q = col->getQOM(iComEl);
  // cout << "Electron. getQOM gives q: "<< q << endl;
  double q = part[iComEl].getQ(pidx);
  // cout << "Electron. getQgives q: "<< q << endl;

  /* Rotates the velocity of the secondary electron */
  secElecVelocity(u, v, w);
  eImpIoniPls[iComEl].create_new_particle(u,v,w,q,x,y,z); 
  // q Here is not QOM -> SHould be +-1 with normalisation in code. 
  // See particles3D.cpp 

  // IoniPls[iComEl].create_new_particle(0.,0.,0.,q,x,y,z);
  // cout << "Completed creating new IoniPls \n";
  // cout << "\n Retrieving info on new IoniPls\n";
  
  /* Save ionospheric ions produced by collisions */
  // Currently ions produced stationary - could be altered.
  // q = col->getQOM(iComIon);
  // cout << "Ion. getQOM gives q: "<< q << endl;
  // q = part[iComIon].getQ(0);
  q = -q;
  // cout << "Ion. getQgives q: "<< q << endl;
  eImpIoniPls[iComIon].create_new_particle(0.,0.,0.,q,x,y,z);
  // Could be given low thermal energy ~0.1 eV
  // COuld/should be done in same way as for photoionization

  // int newidx = eImpIoniPls[iComEl].getNOP();
  // x = eImpIoniPls[iComEl].getX(newidx);
  // y = eImpIoniPls[iComEl].getY(newidx);
  // z = eImpIoniPls[iComEl].getZ(newidx);
  
  // u = eImpIoniPls[iComEl].getU(newidx);
  // v = eImpIoniPls[iComEl].getV(newidx);
  // w = eImpIoniPls[iComEl].getW(newidx);
  // cout << "\n Print properties of saved pls\n ";
  // cout << "x: " << x << " " << y << " " << z << "\n";
  // cout << "v: " << u << " " << v << " " << w << "\n";
  // cout << "Successfully saved Ionized pls.";
}

/* Create Ionized Particles from collisions */
void Collisions::createIonizedParticles(Particles3D *part)
{
  // cout << "\n Creating new particles from Ioni \n";
  for (int i = 0; i < ns; i++)
  {
    // Create new particles in each of the relevant species
    // Only producing cometary electrons/ions
    
    // Producing particles from saved properties in eImpIoniPls
    int nPls = eImpIoniPls[i].getNOP();
    for (int pidx = 0; pidx < nPls; pidx++)
    {
      // Retrieve Properties of secondary pls
      double u = eImpIoniPls[i].getU(pidx);
      double v = eImpIoniPls[i].getV(pidx);
      double w = eImpIoniPls[i].getW(pidx);

      double x = eImpIoniPls[i].getX(pidx);
      double y = eImpIoniPls[i].getY(pidx);
      double z = eImpIoniPls[i].getZ(pidx);
      double q = eImpIoniPls[i].getQ(pidx);
      // cout << "\n Print properties of new pls - Species " << i << " \n ";
      // cout << "x: " << x << " " << y << " " << z << "\n";
      // cout << "v: " << u << " " << v << " " << w << "\n";
    
      // Create secondary particles
      part[i].create_new_particle(u,v,w,q,x,y,z);
      // To add something related to pl ID or tag
      // Want to know they are secondaries, and which parent species
    }
    

  }
  //Clean Ionized particle structures
  resetIoniParticles(); 
}

/* Reset storage of ionized pls. for each timestep*/
void Collisions::resetIoniParticles()
{
  // cout <<"\n Resetting Ioni Pls struct \n";
  for (int i = 0; i < ns; i++)
  {
    // cout << "Species No: " << i << "\n";

    int nPls = eImpIoniPls[i].getNOP();
    // cout << "Species pls: " << nPls <<  "\n";
    for (int pidx = nPls-1; pidx >=0; pidx--) 
    {
      // cout << "pidx: " << pidx << "\n";
      eImpIoniPls[i].delete_particle(pidx);
    }
    // nPls = eImpIoniPls[i].getNOP();
    // cout << "\n No pls in species after deletion: "<< nPls << "\n";
  }
  // cout << "\n All Ioni Pls. deleted\n";
  
}

/* Calculate velocity of secondary pl with random direction 
  Takes velocity of incident electron (post collision as input) */
void Collisions::secElecVelocity(double &u, double &v, double &w)
{
  // two random numbers for theta and phi directions. 
  double rTheta = sample_clopen_u_double(); 
  double rPhi = sample_clopen_u_double();
  // Convert random numbers to angles in radians.
  double theta = acos(2 * rTheta - 1);
  double phi = 8 * atan(1.0) * rPhi; //2*pi* rPhi

  // Calculate velocity magnitude
  double vMag = sqrt(pow(u,2) + pow(v,2) + pow(w,2));
  u = vMag * sin(theta);
  v = u * sin(phi);
  u = u * cos(phi);
  w = vMag * cos(theta);

  // cout << "New pl velocity: " << u << " " << v << " " << w << endl;
}