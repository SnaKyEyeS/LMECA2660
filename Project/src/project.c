#include "poisson.h"
#include "airfoil.h"
#include "cylinder.h"
#include "mesh.h"
#include "stdio.h"
#include "solver.h"
#include "debug.h"



int main(int argc, char *argv[]){


    PetscInitialize(&argc, &argv, 0, 0);

    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(CYLINDER);
    IterateCache *ic = initIterateCache(mesh);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh);

    int debug = 0;

    if (debug) {
        printf("Entered debugging mode\n");
        run_tests();
        exit(1);
    }


    double state = 0.0;         // time
    double dt    = mesh->dt;    // detla-time
    double nu    = mesh->nu;

    double endState;
    switch (mesh->type) {
        case AIRFOIL:
            endState = 15 * mesh->Lc / mesh->Uinf;
            break;

        case CYLINDER:
            endState = 50 * mesh->Lc / mesh->Uinf;
            break;

        default:
            printf("Unknown mapping type.\n");
            exit(1);
    }

    int every_n = 25;
    int max_n = ceil(endState / dt);

    printf("Opening files\n");
    Mesh *meshes[N_MESH] = {mesh->w, mesh->u, mesh->v, mesh->p};

    // Only save W ? //
    int only_w = 1;
    //---------------//

    int N_SAVES = only_w ? 1 : N_MESH;

    FILE *files[N_MESH]  = {
        fopen("../data/mesh_w.txt", "w+"),
        only_w ? NULL : fopen("../data/mesh_u.txt", "w+"),
        only_w ? NULL : fopen("../data/mesh_v.txt", "w+"),
        only_w ? NULL : fopen("../data/mesh_p.txt", "w+")
    };
    
    FILE *diag = fopen("../data/diagnostics.txt", "w+");
    printf("Writing headers and initial states\n");
    for (int i = 0; i < N_SAVES; i++) {
        save_mesh_header(meshes[i], files[i]);
        save_mesh_state(meshes[i], state, files[i]);
    }

    double re, cd, cl, y_plus, max_uv;

    printf("Beginning iterations\n");
    while (state < endState) {
        printf("\tIterate tU/D = %.5f... [%4d/%d]\n", state * mesh->Uinf / mesh->Lc, ic->n+1, max_n);

        iterate(mesh, poisson, ic, state);
        compute_diagnostics(mesh, &cd, &cl, &re, &y_plus, &max_uv, true);

        state += dt;
        if (ic->n % every_n == 0) {
            printf("Saving state.\n");
            for (int i = 0; i < N_SAVES; i++) {
                save_mesh_state(meshes[i], state * mesh->Uinf / mesh->Lc, files[i]);
                fprintf(diag, "%f, %f, %f, %f, %f, %f\n", state * mesh->Uinf / mesh->Lc, cd, cl, re, y_plus, max_uv);
            }
        }
        printf("\n");
    }

    printf("Saving last state !\n");
    for (int i = 0; i < N_SAVES; i++) {
        save_mesh_state(meshes[i], state * mesh->Uinf / mesh->Lc, files[i]);
        fclose(files[i]);
    }

    // Free memory
    free_mac_mesh(mesh);
    free_poisson_solver(poisson);
    freeIterateCache(ic);
    PetscFinalize();
}
