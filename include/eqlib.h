#pragma once

void 
velocity_eq_scaler(
                  double** velocities, 
                  double tau_T, 
                  double dt, 
                  double temperature, 
                  double T_eq,
                  int n_atoms
                  );

void 
pressure_eq_scaler(
                  double** positions, 
                  double alpha_P_cube_root,
                  int n_atoms
                  );

double get_alpha_P_cube_root(
                            double pressure, 
                            double tau_P,
                            double P_eq, 
                            double dt
                            );
