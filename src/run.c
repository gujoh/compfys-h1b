#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include <gsl/gsl_rng.h>
#include <time.h>

void print_positions(double** positions, int n_atoms);
void get_linspace(double* linspace, double start, double end, int n);
void velocity_verlet_one_step(double** positions, double** velocities, double** force,
    double mass, double dt, int n_atoms, double* potential,
    double* kinetic, double* virial, double cell_length);
void task1(void);
void task2(void);
double displacement_magnitude(void);

int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code

    // TASK 1
    //task1();

    // TASK 2
    task2();

    //print_positions(positions, n_atoms);

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
    fclose(file);
}


void task2(void)
{
    FILE* file = fopen("data/task2.csv", "w+");
    int n = 4; 
    int n_atoms = 256; 
    double mass = 0.00279630417;
    double** positions = create_2D_array(n_atoms, 3);
    double** velocities = create_2D_array(n_atoms, 3);
    double** force = create_2D_array(n_atoms, 3);
    double dt;
    printf("dt: ");
    scanf("%lf", &dt);
    double time = 50;
    double timesteps = time / dt;
    double potential = 0;
    double virial = 0;
    double kinetic = 0;
    double lattice_param = 4.03;
    init_fcc(positions, n, lattice_param);
    for (int i = 0; i < n_atoms; i++)
    {
        double displacement = displacement_magnitude();
        printf("%f\n", displacement);
        for (int j = 0; j < 3; j++)
        {
            positions[i][j] += positions[i][j] * displacement;
        }
    }

    calculate(&potential, &virial, force, positions, n * lattice_param, n_atoms);
    for (int t = 0; t < timesteps; t++)
    {
        kinetic = 0;
        velocity_verlet_one_step(positions, velocities, force,
            mass, dt, n_atoms, &potential, &kinetic, &virial, n * lattice_param);
        fprintf(file, "%f,%f\n", potential, kinetic); 
        if (t % 1000 == 0)
        {
            printf("Progress: %.1f%%\n", (t / (float) timesteps) * 100);
        }
    }
    fclose(file);
}

void velocity_verlet_one_step(double** positions, double** velocities, double** force,
    double mass, double dt, int n_atoms, double* potential,
    double* kinetic, double* virial, double cell_length)
{
    for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocities[i][j] += (force[i][j] / mass) * dt / 2;
            positions[i][j] += velocities[i][j] * dt;
        }
    }
    calculate(potential, virial, force, positions, cell_length, n_atoms);
    for (int i = 0; i < n_atoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocities[i][j] += (force[i][j] / mass) * dt / 2;
        }
        double norm = vector_norm(velocities[i], 3);
        *kinetic += 0.5 * mass * norm * norm;
    }
}

double displacement_magnitude(void){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    time_t seed = time(NULL);
    gsl_rng_set(r, seed);
    double u = (gsl_rng_uniform(r) - 0.5) * 0.13;
    return u;
}
