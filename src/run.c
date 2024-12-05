#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include <gsl/gsl_rng.h>
#include <time.h>
#include "eqlib.h"
#include <string.h>

#define K_B  8.617333262e-5 // Boltzmann constant 
#define GPA_TO_BAR 10000 //GPa to bar
#define BAR_TO_EV_A3 0.00000624 // Bar to eV/Å^3

void print_positions(double** positions, int n_atoms);
void get_linspace(double* linspace, double start, double end, int n);
void velocity_verlet_one_step(double** positions, double** velocities, double** force,
    double mass, double dt, int n_atoms, double* potential,
    double* kinetic, double* virial, double cell_length);
double get_msd(double prev_positions[50][256][3], double** positions, int n_atoms, int n_prev_positions);
void task1(void);
void task2(void);
void task3(void);
void task4(void);
gsl_rng* get_rand(void);
void rand_fcc(double** positions, double lattice_param, int n_atoms);

int
run(
    int argc,
    char *argv[]
   )
{

    //task1();

    //task2();

    task3();

    //task4();

    return 0;
}

void print_positions(double** positions, int n_atoms)
{
 for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("%i: %f\t", j, positions[i][j]);
        }
        printf("\n");
    }
}

void get_linspace(double* linspace, double start, double end, int n) 
{
    double increment = (end - start) / n;
    for (int i = 0; i < n; i++)
    {
        linspace[i] = start + i * increment;
    }
}

void get_dt_time_timestep(double* dt, double* time, double* timestep){

    printf("dt: ");
    scanf("%lf", dt);
    printf("max time: ");
    scanf("%lf", time);
    *timestep = *time / *dt;
}

void task1(void)
{
    int n = 4;
    int n_atoms = 256;
    double start = 4;
    double end = 4.1;
    int linspace_len = 1000;
    double* lattice_params = (double*) malloc(sizeof(double) * linspace_len);
    get_linspace(lattice_params, start, end, linspace_len);
    double** positions = create_2D_array(n_atoms, 3);
   
    FILE* file = fopen("data/task1.csv", "w+");
    for (int i = 0; i < linspace_len; i++)
    {
        double lattice_param = lattice_params[i];
        init_fcc(positions, n, lattice_param);
        double e_pot = get_energy_AL(positions, n * lattice_param, n_atoms);
        fprintf(file, "%f,%f\n", pow(lattice_param, 3), e_pot / pow(n, 3));
    }
    free(lattice_params);
    destroy_2D_array(positions);
    fclose(file);
}

void task2(void)
{

    int n = 4; 
    int n_atoms = 256; 
    double mass = 0.00279630417;
    
    double** positions = create_2D_array(n_atoms, 3);
    double** velocities = create_2D_array(n_atoms, 3);
    double** force = create_2D_array(n_atoms, 3);
    
    double dt, time, timesteps;
    get_dt_time_timestep(&dt, &time, &timesteps);

    char buffer[50];
    sprintf(buffer, "data/task2_%.1e.csv", dt);
    FILE* file = fopen(buffer, "w+");

    double potential = 0;
    double virial = 0;
    double lattice_param = 4.046;

    init_fcc(positions, n, lattice_param);
    rand_fcc(positions, lattice_param, n_atoms);

    // Integrating the system.
    calculate(&potential, &virial, force, positions, n * lattice_param, n_atoms);
    for (int t = 0; t < timesteps; t++)
    {
        double kinetic = 0;
        velocity_verlet_one_step(positions, velocities, force,
            mass, dt, n_atoms, &potential, &kinetic, &virial, n * lattice_param);
        
        // T(t) = \frac{2}{3Nk_b}sum\limits_{i=1}^N \frac{p_i^2(t)}{2m_i}
        double temperature = 2.0 / (3.0 * K_B * n_atoms) * kinetic;
        
        // Potential, kinetic, temperature, pressure
        fprintf(file, "%f,%f,%f\n", potential, kinetic, temperature); 
        
        if (t % 50 == 0)
        {
            printf("Progress: %.1f%%\n", (t / (float) timesteps) * 100);
        }
    }
    destroy_2D_array(positions);
    destroy_2D_array(velocities);
    destroy_2D_array(force);
    fclose(file);
}

void task3(void)
{
    int n = 4; 
    int n_atoms = 256; 
    double mass = 0.00279630417;

    int n_prev_positions = 50;
    double** positions = create_2D_array(n_atoms, 3);
    double** velocities = create_2D_array(n_atoms, 3);
    double** force = create_2D_array(n_atoms, 3);
    double prev_positions[n_prev_positions][n_atoms][3];
    memset(prev_positions, 0, sizeof(double) * n_prev_positions * n_atoms * 3);

    double dt, time, timesteps;
    get_dt_time_timestep(&dt, &time, &timesteps);

    char buffer[50];
    
    sprintf(buffer, "data/task3_%.1e.csv", dt);
    FILE* file = fopen(buffer, "w+");
    
    double potential = 0;
    double virial = 0;
    double lattice_param = 4.046; 
    
    double T_eq = 500 + 273.15;
    double P_eq = 1 * BAR_TO_EV_A3;
    double tau_T = 300 * dt;
    double tau_P = 400 * dt;
    int equib_time = 10000;
    
    init_fcc(positions, n, lattice_param); 
    rand_fcc(positions, lattice_param, n_atoms);

    // Integrating the system.
    calculate(&potential, &virial, force, positions, n * lattice_param, n_atoms);
    for (int t = 0; t < timesteps; t++)
    {
        double kinetic = 0;
        velocity_verlet_one_step(positions, velocities, force,
            mass, dt, n_atoms, &potential, &kinetic, &virial, n * lattice_param);
        
        double msd = get_msd(prev_positions, positions, n_atoms, n_prev_positions);
        memcpy(prev_positions[t % n_prev_positions], positions, n_atoms * 3);
        
        // T(t) = \frac{2}{3Nk_b}sum\limits_{i=1}^N \frac{p_i^2(t)}{2m_i}
        double temperature = 2.0 / (3.0 * K_B * n_atoms) * kinetic;
        double volume = n_atoms * pow(lattice_param, 3); 
        double pressure = (n_atoms * K_B * temperature + virial) / volume;

        double alpha_P_cube_root = get_alpha_P_cube_root(pressure, tau_P, P_eq, dt);
        lattice_param *= alpha_P_cube_root;
        if (t < equib_time)
        {
            velocity_eq_scaler(velocities, tau_T, dt, temperature, T_eq, n_atoms);
            pressure_eq_scaler(positions, alpha_P_cube_root, n_atoms);
        }
        // Potential, kinetic, temperature, pressure, lattice constant, positions
        fprintf(file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", potential,
         kinetic, temperature, pressure, lattice_param, 
         positions[60][0], positions[60][1], positions[60][2],
         positions[75][0], positions[75][1], positions[75][2], 
         positions[110][0], positions[110][1], positions[110][2],
         positions[125][0], positions[125][1], positions[125][2],
         msd); 
        
        if (t % 50 == 0)
        {
            printf("Progress: %.1f%%, a = %f, p = %f, temp = %f\n", (t / (float) timesteps) * 100, lattice_param, pressure, temperature);
        }
    }
    destroy_2D_array(positions);
    destroy_2D_array(velocities);
    destroy_2D_array(force);
    fclose(file);
}

void task4(void)
{
    int n = 4; 
    int n_atoms = 256; 
    double mass = 0.00279630417;

    int n_prev_positions = 50;
    double** positions = create_2D_array(n_atoms, 3);
    double** velocities = create_2D_array(n_atoms, 3);
    double** force = create_2D_array(n_atoms, 3);
    double prev_positions[n_prev_positions][n_atoms][3];
    memset(prev_positions, 0, sizeof(double) * n_prev_positions * n_atoms * 3);

    double dt, time, timesteps;
    get_dt_time_timestep(&dt, &time, &timesteps);

    char buffer[50];
    
    sprintf(buffer, "data/task4_%.1e.csv", dt);
    FILE* file = fopen(buffer, "w+");
    
    double potential = 0;
    double virial = 0;
    double lattice_param = 4.046; 
    
    double T_eq = 1000 + 273.15;
    double P_eq = 1 * BAR_TO_EV_A3;
    double tau_T = 300 * dt;
    double tau_P = 400 * dt;
    int equib_time = 10000;
    int equib_time2 = 15000;
    
    init_fcc(positions, n, lattice_param); 
    rand_fcc(positions, lattice_param, n_atoms);

    // Integrating the system.
    calculate(&potential, &virial, force, positions, n * lattice_param, n_atoms);
    for (int t = 0; t < timesteps; t++)
    {
        double kinetic = 0;
        velocity_verlet_one_step(positions, velocities, force,
            mass, dt, n_atoms, &potential, &kinetic, &virial, n * lattice_param);

        double msd = get_msd(prev_positions, positions, n_atoms, n_prev_positions);
        memcpy(prev_positions[t % n_prev_positions], positions, n_atoms * 3);
        
        // T(t) = \frac{2}{3Nk_b}sum\limits_{i=1}^N \frac{p_i^2(t)}{2m_i}
        double temperature = 2.0 / (3.0 * K_B * n_atoms) * kinetic;
        double volume = n_atoms * pow(lattice_param, 3); 
        double pressure = (n_atoms * K_B * temperature + virial) / volume;

        double alpha_P_cube_root = get_alpha_P_cube_root(pressure, tau_P, P_eq, dt);
        lattice_param *= alpha_P_cube_root;
        if (t > equib_time)
        {
            T_eq = 700 + 273.15;
        }
        if (t < equib_time2)
        {
            velocity_eq_scaler(velocities, tau_T, dt, temperature, T_eq, n_atoms);
            pressure_eq_scaler(positions, alpha_P_cube_root, n_atoms);
        }
        // Potential, kinetic, temperature, pressure, lattice constant, positions
        fprintf(file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", potential,
         kinetic, temperature, pressure, lattice_param, 
         positions[60][0], positions[60][1], positions[60][2],
         positions[75][0], positions[75][1], positions[75][2], 
         positions[110][0], positions[110][1], positions[110][2],
         positions[125][0], positions[125][1], positions[125][2],
         msd); 
        
        if (t % 50 == 0)
        {
            printf("Progress: %.1f%%, a = %f, p = %f, temp = %f\n", (t / (float) timesteps) * 100, lattice_param, pressure, temperature);
        }
    }
    destroy_2D_array(positions);
    destroy_2D_array(velocities);
    destroy_2D_array(force);
    fclose(file);
}

// Computes one step of the Velocity verlet integration scheme. 
void velocity_verlet_one_step(double** positions, double** velocities, double** force,
    double mass, double dt, int n_atoms, double* potential,
    double* kinetic, double* virial, double cell_length)
{
    for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocities[i][j] += (force[i][j] / mass) * (dt / 2);
            positions[i][j] += velocities[i][j] * dt;
        }
    }
    calculate(potential, virial, force, positions, cell_length, n_atoms);
    for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocities[i][j] += (force[i][j] / mass) * (dt / 2);
        }
        double norm = vector_norm(velocities[i], 3);
        *kinetic += 0.5 * mass * norm * norm;
    }
}

// Intitializes the random number generator. 
gsl_rng* get_rand(void){
    const gsl_rng_type* T;
    gsl_rng* r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    time_t seed = time(NULL);
    gsl_rng_set(r, seed);
    return r;
}

void rand_fcc(double** positions, double lattice_param, int n_atoms){ 
    gsl_rng* r = get_rand();

    for (int i = 0; i < n_atoms; i++)
    { // Adding randomness to the initial positions. 
        double displacement = (gsl_rng_uniform(r) - 0.5) * 0.13 * lattice_param;
        addition_with_constant(positions[i], positions[i], displacement, 3);
    }
}

double get_msd(double prev_positions[50][256][3], double** positions, int n_atoms, int n_prev_positions)
{
    double msd = 0; 
    for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < n_prev_positions - 1; j++)
        {
            double delta = distance_between_vectors(positions[i], prev_positions[j][i], 3);
            msd += delta * delta;
        }
    }
    return msd / n_atoms / n_prev_positions;
}
