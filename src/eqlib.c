#include "tools.h"
#include <math.h>

void velocity_eq_scaler(double** velocities, double tau_T, double dt, double temperature, double T_eq, int n_atoms){

	double alpha_T = sqrt(1 + (2 * dt * (T_eq - temperature))/(tau_T * temperature));

	for(int i = 0; i < n_atoms; ++i){

		multiplication_with_constant(velocities[i], velocities[i], alpha_T, 3);
	} 
}

void pressure_eq_scaler(double** positions, double alpha_P_cube_root, int n_atoms){
	
	for(int i = 0; i < n_atoms; ++i){

		multiplication_with_constant(positions[i], positions[i], alpha_P_cube_root, 3);
	} 

}

double get_alpha_P_cube_root(double pressure, double tau_P, double P_eq, double dt)
{
	// https://matmake.com/materials-data/aluminum-properties.html
	// Converted from 75.2 GPa to ev/Å^3
	double kappa_T = 1 / 4.69; 
	return pow(1 - kappa_T * (dt/tau_P)*(P_eq - pressure), 1./3.);
}
