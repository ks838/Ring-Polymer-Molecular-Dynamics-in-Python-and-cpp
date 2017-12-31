/*
 @Author:	Kai Song, ks838 _at_ cam.ac.uk

 @Notes :   1. In this edition, we sample from the normal modes 
 			   space directly, without using the coordinate 
 			   transformation too many times (as presented in Ceriotti2010).
            2. In the file 'Eigen_engine' are the necessary parts copied from the 
               famous C++ eigen-3.3.4. 
@Refs   :   1. Craig and David Manolopoulos, J. Chem. Phys. 121, 22  2004
		    2. Ceriotti et al. J. Chem. Phys. 133, 124104  2010 
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <utility>
#include <time.h>

#include "./Eigen_engine/Dense"

#include "consts_rpmd.h"

using namespace std;
using namespace Eigen;

// VectorXd set_pot(VectorXd q){
//	VectorXd V(n_beads);
//	for(int i_bead=0;i_bead<n_beads;++i_bead){
//		V(i_bead) = -a4*pow(q(i_bead),4);  
//	}
//	return V;
//}
clock_t time_start = clock();

VectorXd set_force(VectorXd q){
	VectorXd F(n_beads);
	for(int i_bead=0;i_bead<n_beads;++i_bead){
		//F(i_bead) = -a2*2*pow(q(i_bead),1);  
		F(i_bead) = -a4*4*pow(q(i_bead),3);  

	}
	return F;	
}

MatrixXd trans_matrix(){
	MatrixXd C_matrix(n_beads,n_beads);
	for(int j_bead=0;j_bead<n_beads;j_bead++){
		for(int k_bead=0;k_bead<n_beads;k_bead++){
			if(k_bead == 0){
				C_matrix(j_bead,k_bead) = sqrt(1.0/n_beads);
			}
			else if((k_bead >= 1)&&(k_bead <= (n_beads/2-1))){
				C_matrix(j_bead,k_bead) = sqrt(2.0/n_beads)*cos(2*pi*j_bead*k_bead/n_beads);
			}
			else if(k_bead == n_beads/2){
				C_matrix(j_bead,k_bead) = sqrt(1.0/n_beads)*pow(-1,j_bead);
			}
			else{
				C_matrix(j_bead,k_bead) = sqrt(2.0/n_beads)*sin(2*pi*j_bead*k_bead/n_beads);
			}
		}
	}
	return C_matrix;
}
Matrix2d evol_matrix(double omegak){
	Matrix2d E_matrix(2,2);
	// E_matrix, dt is a global parameter, defined in consts_rpmd
	if(omegak == 0.0){
		E_matrix(0,0) = 1.0;
		E_matrix(0,1) = 0;
		E_matrix(1,0) = dt/mass;
		E_matrix(1,1) = 1.0;
	}
	else{
		E_matrix(0,0) = cos(omegak*dt);
		E_matrix(0,1) = -mass*omegak*sin(omegak*dt);
		E_matrix(1,0) = 1.0/(mass*omegak)*sin(omegak*dt);
		E_matrix(1,1) = cos(omegak*dt);
	}
	return E_matrix;
}

// -------- Langevin Thermostatting------------
pair<VectorXd,VectorXd> C_1_2(VectorXd omegak){
	VectorXd gamma_vec(n_beads),c_1(n_beads),c_2(n_beads);
	gamma_vec = VectorXd::Constant(n_beads,0);
	c_1       = VectorXd::Constant(n_beads,0);
	c_2       = VectorXd::Constant(n_beads,0);

	gamma_vec(0) = 1.0/tau0;
	for(int i_bead=1;i_bead<n_beads; ++i_bead){
		gamma_vec(i_bead) = 2.0*omegak(i_bead);
	}
	for(int i_bead=0;i_bead<n_beads;++i_bead){
		c_1(i_bead) = exp(-dt/2.0*gamma_vec(i_bead));
		c_2(i_bead) = sqrt(1.0-pow(c_1(i_bead),2));
	}
	return make_pair(c_1, c_2);
}


// =========================================================================
//                 ----------EOM_Symplectic--------
// =========================================================================
// @Refs:    
int main(){
	// we just define here
	VectorXd omegak(n_beads); 
	for(int i_bead=0;i_bead < n_beads;++i_bead){
	omegak(i_bead) = 2* omega_N * sin(i_bead*pi/n_beads);
	}
	
	Matrix2d E_matrix(2,2);
	MatrixXd phase_point(n_beads,2);

	pair<VectorXd,VectorXd> c_pair;//c_1(n_beads), c_2(n_beads);
	c_pair = C_1_2(omegak);
	//cout<<c_pair.first<<endl;
	MatrixXd C_matrix(n_beads,n_beads);
	C_matrix=trans_matrix();
	ofstream fout1("xx.dat");	
	// random
	default_random_engine generator_p;
    normal_distribution<double> p_gauss(mu_p,sigma_p);
	default_random_engine generator_langevin;
    normal_distribution<double> xi_gauss(0,1);    
	//default_random_engine generator_uni;
	random_device rd;  //used to obtain a seed for the random number engine
    std::mt19937 gen(rd());//Standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<double> distribution_uni(0,1);    


	VectorXd xx_correlation(nsteps_dynamics);
	xx_correlation = VectorXd::Constant(nsteps_dynamics,0);

	for(int i_sampling=0;i_sampling<n_samplings;++i_sampling){	

	// **********************************************************************
	//         		     	Step 1:  thermalization 
	// **********************************************************************	

		// 1.1 initialize forces, from q_, by set_force(q_)
		VectorXd p_(n_beads);
		p_ = VectorXd::Constant(n_beads,0);
		VectorXd q_(n_beads);
		q_ = VectorXd::Constant(n_beads,0);	

		for(int i_bead=0;i_bead<n_beads;++i_bead){
			//p_(i_bead) = p_gauss(gen);
			p_(i_bead) = p_gauss(rd);
		}	
		//cout<<"p_= "<<p_<<endl;
		//assert(1>2);
		for(int istep=0;istep<nsteps_equil;istep++){
			// -------------- Langevin thermostatting ----------------
			// transform to normal modes
			phase_point.col(0) = C_matrix.transpose()*p_;
			for (int i_bead=0; i_bead<n_beads;++i_bead){
			phase_point(i_bead,0) = c_pair.first(i_bead)*phase_point(i_bead,0)+sqrt(mass/beta_N)* \
								    c_pair.second(i_bead)*xi_gauss(generator_langevin);
			}
			// transform back to plain modes
			p_ = C_matrix*phase_point.col(0);
			p_ += dt*0.5*set_force(q_);
			// transform to normal modes
			// The transpose $a^T$, conjugate $\bar{a}$, and conjugate transpose $a^*$ 
			//  are obtained by the member functions 
			// transpose(), conjugate(), and adjoint(), respectively.
			phase_point.col(0) = C_matrix.transpose()*p_;
			phase_point.col(1) = C_matrix.transpose()*q_;

	// -------------Evolve normal modes------------
			for(int i_bead=0;i_bead<n_beads;++i_bead){
				E_matrix = evol_matrix(omegak(i_bead));
				phase_point(i_bead,0) = E_matrix(0,0)*phase_point(i_bead,0)+ \
										E_matrix(0,1)*phase_point(i_bead,1);
				phase_point(i_bead,1) = E_matrix(1,0)*phase_point(i_bead,0)+ \
										E_matrix(1,1)*phase_point(i_bead,1);						
			}
	// tansform back to plain momenta and coordinates
			p_ = C_matrix*phase_point.col(0);
			q_ = C_matrix*phase_point.col(1);
			p_ += dt/2.0*set_force(q_);
			// -------------- Langevin thermostatting ----------------
			// transform to normal modes
			phase_point.col(0) = C_matrix.transpose()*p_;
			for (int i_bead=0; i_bead<n_beads;++i_bead){
			phase_point(i_bead,0) = c_pair.first(i_bead)*phase_point(i_bead,0)+sqrt(mass/beta_N)* \
								    c_pair.second(i_bead)*xi_gauss(generator_langevin);
			}
			// transform back to plain modes
			p_ = C_matrix*phase_point.col(0);

//	// -------------Andersen: collision------------
//			for(int i_bead=0;i_bead<n_beads;++i_bead){
//				//cout<<distribution_uni(generator_uni)<<endl;
//				if (distribution_uni(gen) <= dt*nu_poisson){
//					p_(i_bead) = p_gauss(generator_p); 
//					//cout<<i_bead<<'\t'<<p_(i_bead)<<endl;
//				}
//			}
		}	

	// --------------END OF THERMALIZATION---------------
		double A_N = 0.0;
		for(int i_bead=0;i_bead<n_beads;++i_bead){
			A_N += q_(i_bead);
		}
		A_N /= n_beads; // A_N(0)

	//cout<<"THERMALIZATION FINISHED ~ ~"<<endl;	

	// ***************************************************************************
	// 				   Step 2:     extracting information
	// ***************************************************************************
		// 2.1 initialize forces, from q_, by set_force(q_)	

		// initial states for the EOMs
		for(int istep=0; istep<nsteps_dynamics;istep++){
			p_ += dt*0.5*set_force(q_);
			// transform to nornal modes
			phase_point.col(0) = C_matrix.transpose()*p_;
			phase_point.col(1) = C_matrix.transpose()*q_;

	//-----------EVOLUTION with normal modes------------
			for(int i_bead=0; i_bead<n_beads;++i_bead){
//				phase_point.row(i_bead) = evol_matrix(omegak(i_bead))* \
											phase_point.row(i_bead);
				E_matrix = evol_matrix(omegak(i_bead));
				phase_point(i_bead,0) = E_matrix(0,0)*phase_point(i_bead,0)+ \
										E_matrix(0,1)*phase_point(i_bead,1);
				phase_point(i_bead,1) = E_matrix(1,0)*phase_point(i_bead,0)+ \
										E_matrix(1,1)*phase_point(i_bead,1);								
			}					
	// tansform back to plain momenta and coordinates
			p_ = C_matrix*phase_point.col(0);
			q_ = C_matrix*phase_point.col(1);

			p_ += dt/2.0*set_force(q_);	

			double B_N = 0.0;
			for(int i_bead =0; i_bead <n_beads;++i_bead){
				B_N += q_(i_bead);
			}
			B_N /= n_beads;
			//print('B_N = ',B_N)
			xx_correlation(istep) += A_N*B_N;
		}
	}// i_sampling	

	for(int istep=0;istep<nsteps_dynamics; istep++){
		fout1<<setw(10)<<dt*istep<<setw(15)<<xx_correlation(istep)/n_samplings<<endl;
	}
	fout1.close();
	printf("Time taken: %.2fs\n", (double)(clock()-time_start)/CLOCKS_PER_SEC);
	return 0;
}



