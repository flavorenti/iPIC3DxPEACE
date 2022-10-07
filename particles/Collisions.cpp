/*******************************************************************************************
  Collisions.h -  Collisional processes for particles.
  -------------------
developers: Peter Stephenson
 ********************************************************************************************/


#include <iostream>
#include <math.h>
#include "Particles3D.h"
#include "Collective.h"
#include "EMfields3D.h"
#include "Collisions.h"
#include "ipicmath.h"
#include "ipicfwd.h"


Collisions::~Collisions()
{

}

Collisions::Collisions(CollectiveIO * col, VirtualTopology3D *vct, Grid * grid, EMfields3D *Emf)
{
// Load simulation parameters
Parameters::init_parameters();
// Load cross sections for simulation
xSec = col->getxSec();

// Species where secondary particles created
iSecElec = col->getiSecElec();
iSecIon = col->getiSecIon();

// No. of collisional processes
nCollProcesses = col->getnCollProcesses();
// Number of ionization collisions
nIoniColls = col->getnIoniColls();
// No of steps where colls skipped
collStepSkip = col->getcollStepSkip();

// Threshold energies
E_th_el = new double[nCollProcesses];
for (int i = 0; i < nCollProcesses; i++){
  E_th_el[i] = col->getE_th_el(i);
}


// Timestep
dt = col->getDt();

// Get Cell size
dx = col->getDx();
dy = col->getDy();
dz = col->getDz();

// No. of species
ns = col->getNs();   

// Get planet offset
PlanetOffset = col->getPlanetOffset();
x_center = col->getx_center();
y_center = col->gety_center();
z_center = col->getz_center();

// Get exosphere properties - Only need while they aren't loaded in Particles3D
R = col->getL_square(); // Mercury/planet radius
Nexo_H  = col->getnSurf(0);   // density of exosphere neutrals at the surface (in nsw units)
hexo_H  = col->gethExo(0); // scale length of exosphere

// Set up ionized Pls object
// Allocation of particles
  eImpIoniPls = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int i = 0; i < ns; i++)
  {
    new(&eImpIoniPls[i]) Particles3D(i,col,vct,grid);
  }

} 

// Apply Collisions from different species. 
void Collisions::Collide(int species, Particles3D *part, CollectiveIO * col, EMfields3D *Emf)
{

  double qom = col->getQOM(species); 

  if (qom<0)// If electrons
  {
    int nPls = part[species].getNOP();
    nCollsLikely = 0.0;
    // Scaling of velocity required when using heavier electrons
    // Not needed for ions

    // Alternative options used for validation
    // elMassRat = 1.0;
    elMassRat = sqrt(qom_eReal / qom); // Electrons
    
    // Alternative options used for validation
    // u0 = col-> getU0(species);
    // v0 = col-> getV0(species);
    // w0 = col-> getW0(species);
    // u0 = 0.0;
    // v0 = 0.0;
    // w0 = 0.0;

    for (int pidx = 0; pidx < nPls; pidx++) 
    {
      // Find partilce position
      xpl = part[species].getX(pidx);
      ypl = part[species].getY(pidx);
      zpl = part[species].getZ(pidx);

      // Find neutral density at particle
      double nNeutr = findneutDensity(species, part);

      // Find Bulk velocity at particle
      findBulkVelocity(species, Emf, u0, v0, w0);
      
      double vScaled = velMagnitude_wScale(species, part, pidx, u0, v0, w0);
      double pColl = ProbColl(nNeutr, vScaled);

      double rColl =  sample_clopen_u_double();// Random number for coll


      if (rColl < pColl) // Particle undergoes collision
      {
        CollideElectron(species, part, pidx, col); // Modify electron velocity
      }
    }
    if  (nCollsLikely>0)
      cout << "No. Pls with large pColl:" << nCollsLikely << endl;
  }
}


// Modififes electron velocities in case of collision. 
// Also saves properties of ionized pls. 
void Collisions::CollideElectron(int species, Particles3D *part, int pidx, CollectiveIO * col)
{
 // Find particle Energy and energy loss for each colliding particle
 // Energy values need to be rescaled into code units. 
  double qom = col->getQOM(species);
  upl = part[species].getU(pidx);
  vpl = part[species].getV(pidx);
  wpl = part[species].getW(pidx);

  Epl = 0.5 * (upl*upl + vpl*vpl + wpl*wpl) / abs(qom); // Pl energy 
  Eth = 0;
  
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
      break;
      
    }
  }

  // Could add continuous cooling processes for v low energies.

  // Reduce particle energies 
  Epl_new = Epl - Eth;
  if (Epl_new<0)
  {
    cout << endl;
    cout << "Epl: " << Epl <<endl;
    cout << "Eth: " << Eth << endl;
    cout << "Epl_new: " << Epl_new << endl;
    eprintf("Particle has energy below zero. ABORT!");
  }
  // Determine new particle velocities
  double vFact = sqrt(Epl_new /Epl);
  
  // Scales electron velocities after energy loss
  upl *= vFact;
  vpl *= vFact;
  wpl *= vFact;

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
double Collisions::ProbColl(double nNeutral, double vMag)
{
  // nNeutral - neutral density of water (?) 
  // Should be same as used for photoionization module

  // xSec - Cross section of collision. Set in input file
  // vScale - Scaling ratio to correct for heavy electrons
  //dprintf("nNeutral is %f",nNeutral);
  double tau = nNeutral * xSec * vMag * dt;
  // double tau = nNeutral * xSec * vMag * dt;
  

  /* If steps skipped then need to scale up collision probability
  by an appropriate amount */
  tau = tau * double(collStepSkip);
  double pColl = 1 - exp(-tau);
  if (pColl > 0.1) nCollsLikely += 1.0;

  return pColl;
}

double Collisions::velMagnitude_wScale(int species, Particles3D *part, int pidx, double u0, 
  double v0, double w0)
{
  // Finds velocity magnitude, including a correction for the increased electron mass.
        
  // Find particle velocity
  upl = part[species].getU(pidx);
  vpl = part[species].getV(pidx);
  wpl = part[species].getW(pidx);
  
  // Calculate scaled up velocities
  double uSC = (upl - u0) * elMassRat + u0;
  double vSC = (vpl - v0) * elMassRat + v0;
  double wSC = (wpl - w0) * elMassRat + w0;
  double vScaled_sq = uSC*uSC + vSC*vSC + wSC*wSC;
  double vScaled = sqrt(vScaled_sq);
  
  return vScaled;
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

  // Information input from primary electrons.

  // Position and velocity of ionizing particle.
  double x = part[species].getX(pidx);
  double y = part[species].getY(pidx);
  double z = part[species].getZ(pidx);
  
  double u = part[species].getU(pidx);
  double v = part[species].getV(pidx);
  double w = part[species].getW(pidx);
  

  // record secondary electron
  // Get weight and sign of charge
  double q = part[species].getQ(pidx);

  /* Rotates the velocity of the secondary electron */
  secElecVelocity(u, v, w);


  // Create electron
  eImpIoniPls[iSecElec].create_new_particle(u,v,w,q,x,y,z); 
  // Create Ion
  q = -q; // Use same weight but opposite charge for ions
  eImpIoniPls[iSecIon].create_new_particle(0.,0.,0.,q,x,y,z);

}

/* Create Ionized Particles from collisions */
void Collisions::createIonizedParticles(Particles3D *part)
{

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

  for (int i = 0; i < ns; i++)
  {
    int nPls = eImpIoniPls[i].getNOP();

    for (int pidx = nPls-1; pidx >=0; pidx--) 
    {
      eImpIoniPls[i].delete_particle(pidx);
    }
  }
  
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
  double vMag = sqrt(u*u + v*v + w*w);
  u = vMag * sin(theta);
  v = u * sin(phi);
  u = u * cos(phi);
  w = vMag * cos(theta);

}


void Collisions::findBulkVelocity(int species, EMfields3D *Emf, double &u0, double &v0, double &w0)
// Finds bulk velocity at the partilce position
{

  // Find cell index closest to pl. 
  int iX = int((xpl/dx) + 0.5);
  int iY = int((ypl/dy) + 0.5);
  int iZ = int((zpl/dz) + 0.5);


  // Find Current and Density
  double Jx = Emf->getJxs(iX, iY, iZ, species);
  double Jy = Emf->getJys(iX, iY, iZ, species);
  double Jz = Emf->getJzs(iX, iY, iZ, species);

  double rho = Emf->getRHOns(iX,iY,iZ,species);

  // Find bulk velocity
  if (rho != 0)
  {
    u0 = Jx / rho;
    v0 = Jy / rho; 
    w0 = Jz / rho;
  }
  else
  {
    u0 = 0.0; v0 = 0.0; w0 = 0.0;
  }

}

double Collisions::findneutDensity(int species, Particles3D *part)
{
  
  // Find distance from planet centre
  double xd = xpl - x_center; // x position relative to planet
  double yd = ypl - y_center; // y position relative to planet
  double zd = zpl - z_center-PlanetOffset; // z position relative to planet
  
  double dist_sq = xd*xd+yd*yd+zd*zd;
  double dist    = sqrt(dist_sq);
  
  // Use neutral density function for exosphere
  double nNeutral = part->neutralDensity(Nexo_H, dist, R, hexo_H);
  
  return nNeutral;
}
