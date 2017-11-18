'''
 @Author:	Kai Song, ks838 _at_ cam.ac.uk

 @Notes :   1. In this edition, we sample from the normal modes 
 			space directly, without using the coordinate 
 			transformation too many times (as presented in Ceriotti2010).

@Refs   :   1. Craig and David Manolopoulos, J. Chem. Phys. 121, 22  2004
		    2. Ceriotti et al. J. Chem. Phys. 133, 124104  2010 
'''
from __future__ import print_function

import numpy as np 
import matplotlib.pyplot as plt 

from consts_rpmd import *
from init_condations import *
from potential import *
from matrix_rpmd import *


# --------------- set the correlation functions :A & B ------------------
'''
 The time correlation functions (TCF) <A(0)B(t)> is a most 
 common/important quantity in quantum dynamics. 
 Presently, the TCF is: <q(0)q(t)>
 You could try, let's say, <q^2(0)q^2(t)>. This is very hard for RMPD 
 to get good results as compared to numerically exact methods, such
 as HEOM.
'''
def A(q):
	#return q*q
	return q
def B(q):
	#return	q*q
	return q
	

# =========================================================================
#                  --------- EOM_Symplectic --------
# =========================================================================
if __name__ == '__main__':
	fout1 = open('qq.dat','w')	

	V_pot = Potential()
	M_rpmd = Matrix_RPMD()
	C_matrix = M_rpmd.trans_matrix()
	C_1, C_2 = M_rpmd.C_1_2()	

	phase_point = np.zeros((n_beads,2))
	qq_correlation = np.zeros(nsteps_dynamics)	

	for i_sampling in range(n_samplings):
	# ****************************************************************
	#         		   	Step 1:  thermalization 
	# ****************************************************************	
		# 1.1 initialize forces, from q_, by V_pot.V_pot.set_force(q_)
		p_ = np.zeros(n_beads)
		q_ = np.zeros(n_beads)
		for i_bead in range(n_beads):
			p_[i_bead] = p_gaussian()	

		for istep in range(nsteps_equil):
			# transform to nornal modes
			phase_point[:,0] = np.dot(p_,C_matrix)
	# ----------Ceriotti: Langevin step-----------
			for i_bead in range(n_beads):
				phase_point[i_bead,0] = C_1[i_bead]*phase_point[i_bead,0] + \
									np.sqrt(mass/beta_N)*C_2[i_bead]*xi_gaussian()
			# transform back to plain modes
			p_ = np.dot(C_matrix, phase_point[:,0])	
			p_ += -dt/2.0*V_pot.set_force(q_)
			# transform to normal modes
			phase_point[:,0] = np.dot(p_, C_matrix)
			phase_point[:,1] = np.dot(q_, C_matrix)
	# -------------Evolve normal modes------------
			for i_bead in range(n_beads):
				phase_point[i_bead][:] = np.dot(M_rpmd.evol_matrix(omegak[i_bead]),
											phase_point[i_bead][:])
	# tansform back to plain momenta and coordinates
			p_ = np.dot(C_matrix,phase_point[:,0])
			q_ = np.dot(C_matrix,phase_point[:,1])
			p_ += -dt/2.0*V_pot.set_force(q_)
	# ----------Ceriotti: Langevin step-----------
			# transorm to normal modes, for PILET
			phase_point[:,0] = np.dot(p_, C_matrix)
			for i_bead in range(n_beads):
				phase_point[i_bead,0] = C_1[i_bead]*phase_point[i_bead,0] + \
									np.sqrt(mass/beta_N)*C_2[i_bead]*xi_gaussian()
			p_ = np.dot(C_matrix, phase_point[:,0])		

	#		# -------------Andersen: collision------------
	#		for i_bead in range(n_beads):
	#			if np.random.uniform(0,1) <= dt* nu_poisson:
	#				p_[i_bead] = p_gaussian() 
	# --------------END OF THERMALIZATION---------------
		A_N = 0.0
		for i_bead in range(n_beads):
			A_N += q_[i_bead]
		A_N /= n_beads # A_N(0)

	# **************************************************************************
	# 				   Step 2:     extracting information
	# **************************************************************************
		# 2.1 initialize forces, from q_, by V_pot.set_force(q_)	

		# initial states for the EOMs
		for istep in range(nsteps_dynamics):
			p_ += -dt*0.5*V_pot.set_force(q_)
	#		# transform to nornal modes
			phase_point[:,0] = np.dot(p_, C_matrix)
			phase_point[:,1] = np.dot(q_, C_matrix)
	#-----------EVOLUTION with normal modes------------
			for i_bead in range(n_beads):
				#omegak = 2* omega_N * np.sin(i_bead*np.pi/n_beads)
				phase_point[i_bead][:] = np.dot(M_rpmd.evol_matrix(omegak[i_bead]),
											phase_point[i_bead][:])
	# tansform back to plain momenta and coordinates
			p_ = np.dot(C_matrix,phase_point[:,0])
			q_ = np.dot(C_matrix,phase_point[:,1])
			p_ += -dt/2.0*V_pot.set_force(q_)

			B_N = 0.0
			for i_bead in range(n_beads):
				B_N += q_[i_bead]
			B_N /= n_beads

			qq_correlation[istep] += A_N*B_N #Z_tot is not needed

	t_list =[0]*nsteps_dynamics
	for istep in range(nsteps_dynamics):
		t_list[istep] = dt*istep
		fout1.write('%14.8f,%18.8f'%(dt*istep,qq_correlation[istep]/n_samplings))
		fout1.write('\n')
	fout1.close()	

	plt.plot(t_list,qq_correlation/n_samplings,linewidth=1.0)
	plt.title("TCF of q")
	plt.xlabel("t (a.u.)")
	plt.ylabel("<q(0)q(t)>")
	plt.show()




