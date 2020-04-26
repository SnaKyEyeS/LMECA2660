#include "poisson.h"
#include "airfoil.h"
#include "cylinder.h"
#include "mesh.h"
#include "stdio.h"
#include "solver.h"
#include "debug.h"

typedef struct {
    int n;          // Current iteration

    // Cache
    double *new_h_x, *new_h_y;
    double *old_h_x, *old_h_y;
    double *grad_P_x, *grad_P_y;
    double *grad_Phi_x, *grad_Phi_y;
    double *nu_lapl_u_x, *nu_lapl_u_y;

} IterateCache;

IterateCache *initIterateCache(MACMesh *mesh) {
    IterateCache *ic = (IterateCache *) malloc(sizeof(IterateCache));

    ic->n = 0;

    // WE do calloc to have, by default, zeros at places we will never change (see B.C.)
    ic->new_h_x = (double *) calloc(sizeof(double), mesh->u->n);
    ic->new_h_y = (double *) calloc(sizeof(double), mesh->v->n);

    ic->old_h_x = (double *) calloc(sizeof(double), mesh->u->n);
    ic->old_h_y = (double *) calloc(sizeof(double), mesh->v->n);

    ic->grad_P_x = (double *) malloc(sizeof(double) * mesh->u->n);
    ic->grad_P_y = (double *) malloc(sizeof(double) * mesh->v->n);

    ic->grad_Phi_x = (double *) malloc(sizeof(double) * mesh->u->n);
    ic->grad_Phi_y = (double *) malloc(sizeof(double) * mesh->v->n);

    ic->nu_lapl_u_x = (double *) malloc(sizeof(double) * mesh->u->n);
    ic->nu_lapl_u_y = (double *) malloc(sizeof(double) * mesh->v->n);

    return ic;
}

void freeIterateCache(IterateCache *ic) {
    free(ic->new_h_x);
    free(ic->new_h_y);

    free(ic->old_h_x);
    free(ic->old_h_y);

    free(ic->grad_P_x);
    free(ic->grad_P_y);

    free(ic->grad_Phi_x);
    free(ic->grad_Phi_y);

    free(ic->nu_lapl_u_x);
    free(ic->nu_lapl_u_y);

    free(ic);
}

/*
 * Set the new value of H to be the old value, so we can over-writte new H.
 * Also increment the variable n by 1.
 */
void icUpdateH(IterateCache *ic) {
    double *temp_x, *temp_y;

    temp_x = ic->old_h_x;
    temp_y = ic->old_h_y;

    ic->old_h_x = ic->new_h_x;
    ic->old_h_y = ic->new_h_y;

    ic->new_h_x = temp_x;
    ic->new_h_y = temp_y;

    ic->n += 1;
}

/*
 * Complete one full iteration, ie.:
 * 1) compute u_star
 * 2) compute phi via poisson solver
 * 3) compute u_n+1
 * 4) compute P^n+1
 */
void iterate(MACMesh *mesh, PoissonData *poisson, IterateCache *ic) {
    double *new_h_x, *new_h_y;
    double *old_h_x, *old_h_y;
    double *u_star, *v_star;
    double *u, *v;
    double *grad_P_x, *grad_P_y;
    double *grad_Phi_x, *grad_Phi_y;
    double *nu_lapl_u_x, *nu_lapl_u_y;
    int ind;

    double dt = mesh->dt;
    double nu = NU;
    // u, v, u*, v* assign

    u = mesh->u->val1;
    v = mesh->v->val1;

    u_star = mesh->u->val2;
    v_star = mesh->v->val2;

    // Assign from the cache

    new_h_x = ic->new_h_x;
    new_h_y = ic->new_h_y;

    old_h_x = ic->old_h_x;
    old_h_y = ic->old_h_y;

    grad_P_x = ic->grad_P_x;
    grad_P_y = ic->grad_P_y;

    grad_Phi_x = ic->grad_Phi_x;
    grad_Phi_y = ic->grad_Phi_y;

    nu_lapl_u_x = ic->nu_lapl_u_x;
    nu_lapl_u_y = ic->nu_lapl_u_y;

    // (1) Compute H, grad_P and nu_lapl
    compute_h(mesh, new_h_x, new_h_y);
    compute_diffusive(mesh, nu_lapl_u_x, nu_lapl_u_y, nu);
    compute_grad(mesh, grad_P_x, grad_P_y, PRESSURE);

    // If first iteration, we use euler explicit
    if (ic->n == 0) {
        // u*
        for (int i = 1; i < mesh->u->n1-1; i++) {
            for (int j = 0; j < mesh->u->n2; j++) {
                ind = i*mesh->u->n2 + j;

                u_star[ind] = u[ind] + dt * (-new_h_x[ind] - grad_P_x[ind] + nu_lapl_u_x[ind]);
            }
        }

        // v*
        for (int i = 0; i < mesh->v->n1; i++) {
            for (int j = 0; j < mesh->v->n2; j++) {
                ind = i*mesh->v->n2 + j;

                v_star[ind] = v[ind] + dt * (-new_h_y[ind] - grad_P_y[ind] + nu_lapl_u_y[ind]);
            }
        }
    }
    else {
        // u*
        for (int i = 1; i < mesh->u->n1-1; i++) {
            for (int j = 0; j < mesh->u->n2; j++) {
                ind = i*mesh->u->n2 + j;

                u_star[ind] = u[ind] + dt * (-0.5 * (3 * new_h_x[ind] - old_h_x[ind]) - grad_P_x[ind] + nu_lapl_u_x[ind]);
            }
        }

        // v*
        for (int i = 0; i < mesh->v->n1; i++) {
            for (int j = 0; j < mesh->v->n2; j++) {
                ind = i*mesh->v->n2 + j;

                v_star[ind] = v[ind] + dt * (-0.5 * (3 * new_h_y[ind] - old_h_y[ind]) - grad_P_y[ind] + nu_lapl_u_y[ind]);
            }
        }
    }

    // Fill u* boundaries

    int ind_inner, ind_outer;

    for (int j = 0; j < mesh->u->n2; j++) {
        ind_inner = j;
        ind_outer = (mesh->u->n1-1)*(mesh->u->n2) + j;

        u_star[ind_inner] = u[ind_inner];
        u_star[ind_outer] = u[ind_outer];
    }

    // (2) Poisson solver

    poisson_solver(poisson, mesh);

    // (3) new speeds

    compute_grad(mesh, grad_Phi_x, grad_Phi_y, PHI);

    // u
    for (int i = 1; i < mesh->u->n1-1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            ind = i*mesh->u->n2 + j;

            u[ind] = u_star[ind] - dt * grad_Phi_x[ind];
        }
    }

    // v
    for (int i = 0; i < mesh->v->n1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;

            v[ind] = v_star[ind] - dt * grad_Phi_y[ind];
        }
    }

    // (4) new P

    for (int i = 0; i < mesh->p->n1; i++) {
        for (int j = 0; j < mesh->p->n2; j++) {
            ind = i*mesh->p->n2 + j;

            mesh->p->val1[ind] += mesh->p->val2[ind];
        }
    }

    icUpdateH(ic);
}

void run_tests() {
    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(AIRFOIL);
    IterateCache *ic = initIterateCache(mesh);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh);

    // Speed tracker
    Tracker *t_u = track_mesh(mesh->u);
    FILE *fp_u = fopen("../data/test_u.txt", "w+");
    save_header(t_u, mesh->u->x, mesh->u->y, fp_u);

    Tracker *t_v = track_mesh(mesh->v);
    FILE *fp_v = fopen("../data/test_v.txt", "w+");
    save_header(t_v, mesh->v->x, mesh->v->y, fp_v);
    // Pressure tracker
    Tracker *t_p = track_mesh(mesh->p);
    FILE *fp_p = fopen("../data/test_p.txt", "w+");
    save_header(t_p, mesh->p->x, mesh->p->y, fp_p);

    // Set speed & pressure
    set_speed_pressure(mesh, solenoidal_speed);

    // Save and close speeds
    save_tracking_state(t_u, 0.0, fp_u);
    fclose(fp_u);

    save_tracking_state(t_v, 0.0, fp_v);
    fclose(fp_v);

    // Save and close pressure
    save_tracking_state(t_p, 0.0, fp_p);
    fclose(fp_p);

    // H tracker, content will be save into val2 (!) of u and v
    Tracker *t_h_x = track_mesh(mesh->u);
    FILE *fp_h_x = fopen("../data/test_h_x.txt", "w+");
    save_header(t_h_x, mesh->u->x, mesh->u->y, fp_h_x);

    Tracker *t_h_y = track_mesh(mesh->v);
    FILE *fp_h_y = fopen("../data/test_h_y.txt", "w+");
    save_header(t_h_y, mesh->v->x, mesh->v->y, fp_h_y);

    compute_h(mesh, mesh->u->val2, mesh->v->val2);

    // Save and close h
    save_tracking_state(t_h_x, 0.0, fp_h_x);
    fclose(fp_h_x);

    save_tracking_state(t_h_y, 0.0, fp_h_y);
    fclose(fp_h_y);

    // Laplacian tracker, content will be save into val2 (!) of u and v
    Tracker *t_lapl_x = track_mesh(mesh->u);
    FILE *fp_lapl_x = fopen("../data/test_lapl_x.txt", "w+");
    save_header(t_lapl_x, mesh->u->x, mesh->u->y, fp_lapl_x);

    Tracker *t_lapl_y = track_mesh(mesh->v);
    FILE *fp_lapl_y = fopen("../data/test_lapl_y.txt", "w+");
    save_header(t_lapl_y, mesh->v->x, mesh->v->y, fp_lapl_y);

    compute_diffusive(mesh, mesh->u->val2, mesh->v->val2, 1.0);

    // Save and close lapl
    save_tracking_state(t_lapl_x, 0.0, fp_lapl_x);
    fclose(fp_lapl_x);

    save_tracking_state(t_lapl_y, 0.0, fp_lapl_y);
    fclose(fp_lapl_y);

    // Grap P tracker, content will be save into val2 (!) of u and v
    Tracker *t_grad_p_x = track_mesh(mesh->u);
    FILE *fp_grad_p_x = fopen("../data/test_grad_p_x.txt", "w+");
    save_header(t_grad_p_x, mesh->u->x, mesh->u->y, fp_grad_p_x);

    Tracker *t_grad_p_y = track_mesh(mesh->v);
    FILE *fp_grad_p_y = fopen("../data/test_grad_p_y.txt", "w+");
    save_header(t_grad_p_y, mesh->v->x, mesh->v->y, fp_grad_p_y);

    compute_grad(mesh, mesh->u->val2, mesh->v->val2, PRESSURE);

    // Save and close grap p
    save_tracking_state(t_grad_p_x, 0.0, fp_grad_p_x);
    fclose(fp_grad_p_x);

    save_tracking_state(t_grad_p_y, 0.0, fp_grad_p_y);
    fclose(fp_grad_p_y);

}

int main(int argc, char *argv[]){


    PetscInitialize(&argc, &argv, 0, 0);

    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(CYLINDER);
    IterateCache *ic = initIterateCache(mesh);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh);

    //run_tests();


    double state = 0.0;         // time
    double dt    = mesh->dt;    // detla-time
    double nu    = NU;

    double endState = 1;

    int every_n = 100;
    int max_n = 2000;

    printf("Opening files\n");
    Mesh *meshes[N_MESH] = {mesh->w, mesh->u, mesh->v, mesh->p};
    FILE *files[N_MESH]  = {
        fopen("../data/mesh_w.txt", "w+"),
        fopen("../data/mesh_u.txt", "w+"),
        fopen("../data/mesh_v.txt", "w+"),
        fopen("../data/mesh_p.txt", "w+")
    };

    printf("Writing headers and initial states\n");
    for (int i = 0; i < N_MESH; i++) {
        save_mesh_header(meshes[i], files[i]);
        save_mesh_state(meshes[i], state, files[i]);
    }

    printf("Beginning iterations\n");
    while (state < endState) {
        printf("\tIterate t=%.5fs... [%4d/%d]\n", state, ic->n+1, max_n);
        iterate(mesh, poisson, ic);

        state += dt;

        if (ic->n % every_n == 0) {
            printf("Saving state.\n");
            for (int i = 0; i < N_MESH; i++) {
                save_mesh_state(meshes[i], state, files[i]);
            }

        }

        if (ic->n == max_n) {
            break;
        }
    }

    for (int i = 0; i < N_MESH; i++) {
        fclose(files[i]);
    }

    // Free memory
    free_mac_mesh(mesh);
    free_poisson_solver(poisson);
    freeIterateCache(ic);
    PetscFinalize();

}
