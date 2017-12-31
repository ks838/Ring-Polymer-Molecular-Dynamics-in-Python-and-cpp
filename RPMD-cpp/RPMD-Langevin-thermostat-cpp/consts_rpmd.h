#ifndef RPMD_CPP_CONSTs_RPMD_H_
#define RPMD_CPP_CONSTs_RPMD_H_
#endif //RPMD_CPP_CONSTs_RPMD_H_

#include <complex>
#include "./Eigen_engine/Dense"

using namespace std;
using namespace Eigen;

// Physical constants
const double bolz   = 1.3806503e-23;
const double au2amu = 5.4857990e-4;
const double au2cm  = 219474.63068;
const double au2fs  = 2.418884254e-2;
const double pi     = 3.1415926535897932;
const complex<double> eye = 1i;

// parameters for the system
const double temperature = 300; // Kelvin
const double beta = 1.0; // a.u.
const double mass = 1.0;
                                
// for the potentials: through changing a1-a4, 
// we could get different double-well potential forms
const double a1 = 0;
const double a2 = 0.5;
const double a3 = 0.1;
const double a4 = 0.25;      

// parameters for propagation
const double dt = 0.01;
const int nsteps_equil = 1950;
const int nsteps_dynamics = 1600;

// for the n beads
const int n_beads = 4; // an even number
const double omega_N = n_beads/beta; // hbar = 1
const double beta_N  = beta/n_beads;
//---omegak would be used many times,

//------for Andersen thermostatting-----
const double nu_poisson = 0.07;
//------for Ceriotti thermostatting-----
const double au0 = 0.7; // an input parameter for tuning the efficiency


// for samplings
const double gamma_adiabatic = 1.5; // the adiabatic parameter
const double omega_aidabatic = gamma_adiabatic/beta_N;

// trajectories
const int n_samplings = 122;

const double mu_p = 0.0;
const double sigma_p = sqrt(mass*n_beads/beta);
const double mu_q = 0.0;

// --------- Langevin thermostatting -------
// @Refs:  Ceriotti 2010
const double tau0 = 0.7; // an input parameter for tuning the efficiency
// ---------- Langevin dissipation ---------
// @Refs:  Mark Tuckerman
const double gamma_LD = 0.8;

