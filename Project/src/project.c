#include "poisson.h"
#include "airfoil.h"
#include "cylinder.h"
#include "mesh.h"
#include "stdio.h"


int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);

    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(CYLINDER);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh);


    // Compute Poisson Solution
    poisson_solver(poisson, mesh);
    save_mesh(mesh->p);

    // Free memory
    free_mac_mesh(mesh);
    free_poisson_solver(poisson);
    PetscFinalize();

}
