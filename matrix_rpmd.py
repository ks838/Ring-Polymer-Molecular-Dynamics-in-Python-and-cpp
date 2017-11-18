'''
 @Author:	Kai Song, ks838 _at_ cam.ac.uk
 @Notes :   This part gives the matrices that would be used in our RPMD calculations.
'''
import numpy as np 
from consts_rpmd import *

class Matrix_RPMD:
	# transform to normal modes
	#Refs:   Ceriotti et al. J. Chem. Phys. 133, 124104  2010 
	def trans_matrix(self):
		C_matrix = np.zeros((n_beads,n_beads))
		for j_bead in range(n_beads):
			for k_bead in range(n_beads):
				if k_bead == 0:
					C_matrix[j_bead][k_bead] = np.sqrt(1/n_beads)
				elif k_bead >= 1 and k_bead <= (n_beads/2-1):
					C_matrix[j_bead][k_bead] = np.sqrt(2/n_beads)*np.cos(2*np.pi*j_bead*k_bead/n_beads)
				elif k_bead == n_beads/2:
					C_matrix[j_bead][k_bead] = np.sqrt(1/n_beads)*np.power(-1,j_bead)
				else:
					C_matrix[j_bead][k_bead] = np.sqrt(2/n_beads)*np.sin(2*np.pi*j_bead*k_bead/n_beads)
		return C_matrix
	# the classical evolution matrix (for normal modes)
	def evol_matrix(self,omegak):
		# dt is a global parameter, defined in consts_rpmd
		if omegak == 0:
			E_matrix[0][0] = 1.0
			E_matrix[0][1] = 0
			E_matrix[1][0] = dt/mass
			E_matrix[1][1] = 1.0
		else:
			E_matrix[0][0] = np.cos(omegak*dt)
			E_matrix[0][1] = -mass*omegak*np.sin(omegak*dt)
			E_matrix[1][0] = 1.0/(mass*omegak)*np.sin(omegak*dt)
			E_matrix[1][1] = np.cos(omegak*dt)
		return E_matrix	
	
	#------- thermostatting ------------
	# Langevin thermostatting is used here. The other commonly-used one is Andersen thermostatting
	# Refs:   Ceriotti et al. J. Chem. Phys. 133, 124104  2010 
	# We denote them as c_1, c_2 just as the paper did.
	def C_1_2(self):
		gamma_vec = np.zeros(n_beads)
		c_1 = np.zeros(n_beads)
		c_2 = np.zeros(n_beads)
		gamma_vec[0] = 1.0/tau0
		for i_bead in range(1,n_beads):
			gamma_vec[i_bead] = 2*omegak[i_bead]
		for i_bead in range(n_beads):
			c_1[i_bead] = np.exp(-dt/2*gamma_vec[i_bead])
			c_2[i_bead] = np.sqrt(1.0-c_1[i_bead]**2)
		return c_1, c_2
