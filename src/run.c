#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

void print_positions(double** positions, int n_atoms);
void get_linspace(double* linspace, double start, double end, int n);

int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
    int n = 4;
    int n_atoms = 256;
    double start = 3.9;
    double end = 4.5;
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
        fprintf(file, "%f,%f\n", lattice_param, e_pot);
    }
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
