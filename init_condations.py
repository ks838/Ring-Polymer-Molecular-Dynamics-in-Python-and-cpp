'''
 @Author:	Kai Song, ks838 _at_ cam.ac.uk

 '''
import numpy as np 
from consts_rpmd import *

#-------------------------Initial Conditions-------------------------
# set the initial states for momenta and positions
# sampled from Gaussian distribution
def p_gaussian():
	#mu_p and simga_p have been defined outside
	p_gaussian = np.random.normal(mu_p,sigma_p)
	return p_gaussian

# set the initial states for positions,
# also Gaussian distribution, as we work in the normal modes
def q_gaussian(i_bead):
#  For the positions, it is helpful to think of the ring polymer centroid 
#  as the “classical” coordinate which can be initialized 

	#omegak = 2* omega_N * np.sin(i_bead*np.pi/n_beads)
	sigma_q = np.sqrt(mass*n_beads/beta)/omegak[i_bead]
	q_gaussian = np.random.normal(mu_q,sigma_q)
	return q_gaussian

# Gaussian-distributed random variables. 
# THe "quality" of these random varialbes do make sense sometimes.
# If you want to generate them by yourself, just use Box-Muller transformation
def xi_gaussian():
	xi_gaussian = np.random.normal(0,1)
	return xi_gaussian