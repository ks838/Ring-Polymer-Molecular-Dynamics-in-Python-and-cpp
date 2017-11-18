'''
 @Author:	Kai Song, ks838 _at_ cam.ac.uk
 @Notes :   This part gives the constants and parameters.
'''
import numpy as np 

# Physical constants
bolz = 1.3806503e-23
au2amu = 5.4857990e-4
au2cm = 219474.63068
au2fs = 2.418884254e-2
eye = complex(0.0,1.0)

# parameters for the system
#the inverse temperature 
beta = 1.0 # a.u.

mass = 1.0
                                
# for the potentials
a1 = 0; a2 = 0.5
a3 = 0.0; a4 = 0.25      

# ------ params for propagation ------
dt = 0.01 # time step
# steps for the equilibrating part
nsteps_equil = 2000
# steps for the dynamics
nsteps_dynamics = 1000

# -------- for the n beads ----------
# for simple potential forms (e.g., a double-well form), n_beads <10 are 
# engough. And, for harmonic form, n_beads = 1 is engough
n_beads = 2 # should be an even number in our settings
omega_N = n_beads/beta # we have used hbar = 1
beta_N  = beta/n_beads

#omegak would be used for many times,
# we just define here 
omegak = np.zeros(n_beads)
for i_bead in range(n_beads):
	omegak[i_bead] = 2* omega_N * np.sin(i_bead*np.pi/n_beads)

E_matrix = np.zeros((2,2))	

#------ parameter for Ceriotti thermostatting-----
tau0 = 0.7 # an input parameter for tuning the efficiency
# for samplings
gamma_adiabatic = 1.5 # the adiabatic parameter
omega_aidabatic = gamma_adiabatic/beta_N

# The number of samplings (from the thermostatting).
# Typically, we need ~10^4 to get converged results.
# We started using a small number for testing
n_samplings = 1000

mu_p = 0.0
sigma_p = np.sqrt(mass*n_beads/beta)
mu_q = 0.0


# set the mass vector for the beads
m_ = np.zeros(n_beads)
m_[0] = mass

for i_bead in range(1,n_beads):
	omega_i = 2.0/beta_N*np.sin(i_bead*np.pi/n_beads)
	m_[i_bead] = mass * (omega_i/omega_aidabatic)**2

