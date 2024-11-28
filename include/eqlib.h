#pragma once

void 
velocity_eq_scaler(
                  double** velocities, 
                  double tau_T, 
                  double dt, 
                  double temperature, 
                  int n_atoms
                  );

void 
pressure_eq_scaler(
                  double** positions, 
                  double pressure, 
                  double tau_P, 
                  double P_eq, 
                  int n_atoms
                  );