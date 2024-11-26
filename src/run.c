#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

void print_positions(double** positions, int n_atoms);

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
    double lattice_param = 0.5;
    double** positions = create_2D_array(n_atoms, 3);
    init_fcc(positions, n, lattice_param);
   
    print_positions(positions, n_atoms);

    double e_pot = get_energy_AL(positions, n, n_atoms);
    printf("%f", e_pot);
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