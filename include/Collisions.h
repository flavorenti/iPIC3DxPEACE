/*******************************************************************************************
  Collisions.h -  Collisional processes for particles.
  -------------------
developers: Peter Stephenson
 ********************************************************************************************/
#ifndef _Collisions_
#define _Collisions_
#include "ipicfwd.h"
/*! @file 
 *  Class to manage collisional processes
    Developers: P. Stephenson (Jul 2022)
 */


class Collisions {
  public:
       /*! @brief Constructor for collision object
       */
      ~Collisions();

      /*! @brief Constructor for collision object
      *
      * Creates the collision object and retrieves parameters from the Collective
      * Also creates object where secondary ions and electrons are stored.
      *
      * @param[in] col Collective with simulation parameters
      * @param[in] vct Topology of the simulation
      * @param[in] grid Grid object
      * @param[in] EMf Fields object 
      */
      Collisions(CollectiveIO * col, VirtualTopology3D *vct, Grid * grid, EMfields3D *Emf);

      /*! @brief Determines if collisions occur and implements collisions.
      *
      * Cycles through particles of each species, testing whether they undergo
      * collisions or not. Energy degradation and recording of secondary particles 
      * from impact ionization occur within.
      * 
      * @param[in] species index of species
      * @param[in] part particles object
      * @param[in] col Collective with simulation parameters
      * @param[in] EMf Fields object 
      */
      void Collide(int species, Particles3D *part, CollectiveIO * col, EMfields3D *Emf);

      /*! @brief Carries out collision for electrons
      *
      * Cycles through particles of each species, testing whether they undergo
      * collisions or not. Energy degradation and recording of secondary particles 
      * from impact ionization occur within.
      * 
      * @param[in] species index of species
      * @param[in] part particles object
      * @param[in] pidx index of colliding particle
      * @param[in] col Collective with simulation parameters
      * @param[in] EMf Fields object 
      */
      void CollideElectron(int species, Particles3D *part, int pidx, CollectiveIO *col);
          // * @param[in] vMag magnitude of colliding electron velocity
      // void CollideElectron(int species, Particles3D *part, int pidx, double vMag,
      //   CollectiveIO *col);

      /*! @brief Calculates the probability of a collision occuring
      *
      * Calculates probability that a particle undergoes a collision in a
      * single step. For electrons, uses a rescaled velocity to correct
      * for the reduced electron-ion mass ratio
      * 
      * @param[in] vMag magnitude of the rescaled velocity
      */
      double ProbColl(double vMag);

      /*! @brief Calculates velocity magnitude
      *
      * Calculates magnitude of the particle velocity, with no scaling
      * to correct for the electron mass.
      * 
      * @param[in] species index of species
      * @param[in] part particles object
      * @param[in] pidx index of colliding particle
      */
      double velMagnitude(int species, Particles3D *part, int pidx);

      /*! @brief Calculates velocity magnitude with correction for electron mass
      *
      * Calculates magnitude of the particle velocity, with a correction
      * for the mass of the electron. If electron, scales up the thermal velocity to correct 
      * for the increased electron mass. 
      * 
      * @param[in] species index of species
      * @param[in] part particles object
      * @param[in] pidx index of colliding particle
      * @param[in] u0 bulk velocity at pl along x
      * @param[in] v0 bulk velocity at pl along y
      * @param[in] w0 bulk velocity at pl along z
      */
      double velMagnitude_wScale(int species, Particles3D *part, int pidx, 
        double u0, double v0, double w0);
      // double velMagnitude_wScale(int species, Particles3D *part, int pidx, 
      //   double u0, double v0, double w0, double u0Mag2, double vScale);

      /*! @brief Calculates velocity magnitude with correction for electron mass
      *
      * NOT DONE. SHOULD BE SAME FUNCTION USED TO GENERATE EXOSPHERE.
      */
      double neutralDensity(Particles3D part, int pidx);

      /*! @brief Records properties of the secondary particles.
      *
      * Records properties of the secondary particles that are produced in 
      * impact-ionization collisions. 
      * 
      * @param[in] species index of species
      * @param[in] part particles object
      * @param[in] pidx index of colliding particle
      * @param[in] col Collective with simulation parameters
      */
      void recordIonizedParticles(int species, Particles3D *part, int pidx, CollectiveIO * col);

      /*! @brief Creates particles generated from impact ionization.
      *
      * Creates secondary particles in the exospheric species, using the saved 
      * particle properties in eImpIoniPls.
      * 
      * @param[inout] part particles object
      */
      void createIonizedParticles(Particles3D *part);

      /*! @brief Creates particles generated from impact ionization.
      *
      * Removes all recorded secondary particle properties after the 
      * particles have been produced. 
      * 
      */
      void resetIoniParticles();

      /*! @brief Finds velocity of secondary electron in random direction
      *
      * Calculates velocity of secondary electron with random angles
      * theta and phi. Uniform angular distribution.
      * 
      * @param[inout] u electron velocity in x-direction
      * @param[inout] v electron velocity in y-direction
      * @param[inout] w electron velocity in z-direction
      */
      void secElecVelocity(double &u, double &v, double &w);
      double velocityScale(Collective *col, int species, double qom);

      /*! @brief Finds bulk velocity at position of electron 
      *
      * Calculates velocity of secondary electron with random angles
      * theta and phi. Uniform angular distribution.
      * 
      * @param[in] species index of species
      * @param[in] EMf Fields object 
      * @param[in] part particles object
      * @param[in] col Collective with simulation parameters
      * @param[in] pidx index of particle
      * @param[out] u0 bulk velocity in x at particle
      * @param[out] v0 bulk velocity in y at particle
      * @param[out] w0 bulk velocity in z at particle
      */
      void findBulkVelocity(int species, EMfields3D *Emf, Particles3D *part, Collective *col,
        int pidx, double &u0, double &v0, double &w0);
 
  private:
    Particles3D *eImpIoniPls; // Particles generated by e-impact ionization
    /*! Collision cross Section */
    double xSec;// Collision cross section in code units. 
    int nCollProcesses; // Number of collisional processes
    double *E_th_el; // Threshold energies for electron collisions.
    int nIoniColls; // No. of ionization collisions
    int collStepSkip;
    bool isIoni;
    /*! Threshold energy for each process*/
 
    double nNeutral = 1e7; //Hard coded neutral density for validation
    double qom_eReal = -1836; // Real charge per mass for electrons (if qom = 1 for protons)
    // Needs to be overuled by function (from photoionzation module)
    // double mPl = 2e4;

    // Used to keep track of how many electrons have 'large' coll probabilities
    double nCollsLikely;
    
    int ns;
    double dt;
    // // Particle Velocity 
    double upl;
    double vpl;
    double wpl;
    // Particle Position
    double xpl;
    double ypl;
    double zpl;
    // Particle energy
    double Epl;
    double Eth;
    double Epl_new;

    // Bulk Velocity
    double u0;
    double v0;
    double w0;

    // Ratio between electron mass in PiC sim. and real world
    // sqrt(me_PiC / me_Real)
    double elMassRat;

    // Species indexes for secondary ions and electrons
    int iSecElec;
    int iSecIon;

  protected:
    const Collective * col;   
    const VirtualTopology3D * vct;
    const Grid * grid;


    
};
#endif