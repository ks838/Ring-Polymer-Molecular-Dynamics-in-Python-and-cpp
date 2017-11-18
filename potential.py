'''
 @Author:	Kai Song, ks838 _at_ cam.ac.uk

'''
from consts_rpmd import *

class Potential:
	# Here, for first-stage testing, we use the simplest (also the most classical) one: a harmonic form.
	# In practice, this can be rather complex, e.g., simulated from Gaussian or molecular dynamics.
	def set_pot(self,q):
		V = a2*q**2  
		return V	
	# We always try to get the analytic derivative of potential V(q) 
	def set_force(self,q):
		F = a2*2*q  
		return F
