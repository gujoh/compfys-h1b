#include "tools.h"
#include <math.h>

void velocity_eq_scaler(double** velocities, double tau_T, double dt, double temperature, int n_atoms){

	double alpha_T = sqrt(1 + (2 *dt * (T_eq - temperature))/(tau_T * temperature));

	for(int i = 0; i < n_atoms; ++i){

		multiplication_with_constant(velocities[i], velocities[i], alpha_T, 3)
	} 
}

void pressure_eq_scaler(double** positions, double pressure, double tau_P, double P_eq, int n_atoms){

	double kappa_T = 1. / (76e3); // (bulk modulus)^{-1}  
	double alpha_P = pow(1 - kappa_T * (dt/tau_P)*(P_eq - pressure), 1/3);
	
	for(int i = 0; i < n_atoms; ++i){

		multiplication_with_constant(positions[i], positions[i], alpha_P, n_atoms);
	} 

}