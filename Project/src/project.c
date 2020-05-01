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
        printf("First iteration with Euleur explicit.\n");
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

    // Compute Re_w,max
    double x, y, r, h;
    double Re_w, Re_w_max = -1;
    int ind_max = 0;
    for (int ind = 0; ind < mesh->w->n; ind++) {
        x = mesh->w->x[ind];
        y = mesh->w->y[ind];
        r = hypot(x, y);

        // r <= 12 * D (D = 1)
        if (r <= 12.0/50.0) {
            h = fmax(mesh->w->h1[ind]*mesh->w->d1, mesh->w->h2[ind]*mesh->w->d2);
            Re_w = fabs(mesh->w->val1[ind]) * h*h / NU;
            if (Re_w > Re_w_max) {
                ind_max = ind;
                Re_w_max = Re_w;
            }
        }
    }

    x = mesh->w->x[ind_max];
    y = mesh->w->y[ind_max];
    r = hypot(x, y);
    printf("\tMesh Reynolds = %f at r = %f\n", Re_w_max, r);

}

void run_tests() {
    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(CYLINDER);
    IterateCache *ic = initIterateCache(mesh);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh);

    /////////////     1     /////////////
    printf("1) Tracking initialisation of u, v and p\n");

    // Speed tracker
    Tracker *t_u = track_mesh(mesh->u, mesh->dt);
    FILE *fp_u = fopen("../data/test_u.txt", "w+");
    save_header(t_u, mesh->u->x, mesh->u->y, fp_u);

    Tracker *t_v = track_mesh(mesh->v, mesh->dt);
    FILE *fp_v = fopen("../data/test_v.txt", "w+");
    save_header(t_v, mesh->v->x, mesh->v->y, fp_v);
    // Pressure tracker
    Tracker *t_p = track_mesh(mesh->p, mesh->dt);
    FILE *fp_p = fopen("../data/test_p.txt", "w+");
    save_header(t_p, mesh->p->x, mesh->p->y, fp_p);

    // UNCOMMENT TO : Set speed & pressure to other function
    // set_speed_pressure(mesh, solenoidal_speed);

    // Save and close speeds
    save_tracking_state(t_u, 0.0, fp_u);
    fclose(fp_u);
    free_Tracker(t_u);

    save_tracking_state(t_v, 0.0, fp_v);
    fclose(fp_v);
    free_Tracker(t_v);

    // Save and close pressure
    save_tracking_state(t_p, 0.0, fp_p);
    fclose(fp_p);
    free_Tracker(t_p);

    printf("1) u, v and p trackers are saved and closed!\n");

    /////////////     2     /////////////
    printf("2) Tracking numerical methods individually\n");
    printf("\t2.1) Tracking H\n");
    // H tracker, content will be save into val2 (!) of u and v
    Tracker *t_h_x = track_mesh(mesh->u, mesh->dt);
    FILE *fp_h_x = fopen("../data/test_h_x.txt", "w+");
    save_header(t_h_x, mesh->u->x, mesh->u->y, fp_h_x);

    Tracker *t_h_y = track_mesh(mesh->v, mesh->dt);
    FILE *fp_h_y = fopen("../data/test_h_y.txt", "w+");
    save_header(t_h_y, mesh->v->x, mesh->v->y, fp_h_y);

    compute_h(mesh, mesh->u->val2, mesh->v->val2);

    // Save and close h
    save_tracking_state(t_h_x, 0.0, fp_h_x);
    fclose(fp_h_x);
    free_Tracker(t_h_x);

    save_tracking_state(t_h_y, 0.0, fp_h_y);
    fclose(fp_h_y);
    free_Tracker(t_h_y);

    printf("\t2.1) H trackers are saved and closed!\n");

    printf("\t2.2) Tracking Laplacian and omega\n");
    // Laplacian tracker, content will be save into val2 (!) of u and v
    Tracker *t_lapl_x = track_mesh(mesh->u, mesh->dt);
    FILE *fp_lapl_x = fopen("../data/test_lapl_x.txt", "w+");
    save_header(t_lapl_x, mesh->u->x, mesh->u->y, fp_lapl_x);

    Tracker *t_lapl_y = track_mesh(mesh->v, mesh->dt);
    FILE *fp_lapl_y = fopen("../data/test_lapl_y.txt", "w+");
    save_header(t_lapl_y, mesh->v->x, mesh->v->y, fp_lapl_y);

    Tracker *t_w = track_mesh(mesh->w, mesh->dt);
    FILE *fp_w = fopen("../data/test_w.txt", "w+");
    save_header(t_w, mesh->w->x, mesh->w->y, fp_w);

    compute_diffusive(mesh, mesh->u->val2, mesh->v->val2, 1.0);

    // Save and close lapl
    save_tracking_state(t_lapl_x, 0.0, fp_lapl_x);
    fclose(fp_lapl_x);
    free_Tracker(t_lapl_x);

    save_tracking_state(t_lapl_y, 0.0, fp_lapl_y);
    fclose(fp_lapl_y);
    free_Tracker(t_lapl_y);

    save_tracking_state(t_w, 0.0, fp_w);
    fclose(fp_w);
    free_Tracker(t_w);

    printf("\t2.2) Laplacian and omega trackers are saved and closed!\n");

    printf("\t2.3) Tracking grad P\n");
    // Grap P tracker, content will be save into val2 (!) of u and v
    Tracker *t_grad_p_x = track_mesh(mesh->u, mesh->dt);
    FILE *fp_grad_p_x = fopen("../data/test_grad_p_x.txt", "w+");
    save_header(t_grad_p_x, mesh->u->x, mesh->u->y, fp_grad_p_x);

    Tracker *t_grad_p_y = track_mesh(mesh->v, mesh->dt);
    FILE *fp_grad_p_y = fopen("../data/test_grad_p_y.txt", "w+");
    save_header(t_grad_p_y, mesh->v->x, mesh->v->y, fp_grad_p_y);

    compute_grad(mesh, mesh->u->val2, mesh->v->val2, PRESSURE);

    // Save and close grap p
    save_tracking_state(t_grad_p_x, 0.0, fp_grad_p_x);
    fclose(fp_grad_p_x);
    free_Tracker(t_grad_p_x);

    save_tracking_state(t_grad_p_y, 0.0, fp_grad_p_y); // Higher error on this term ?
    fclose(fp_grad_p_y);
    free_Tracker(t_grad_p_y);

    printf("\t2.3) Grad P trackers are saved and closed!\n");

    /////////////     3     /////////////
    // Debuging iterations

    printf("3) Tracking first iteration\n");

    printf("dt = %.10f\n", mesh->dt);

    // 1) u*, v* tracker

    printf("\t3.1) Tracking u*, v*\n");
    Tracker *t_u_star = track_mesh(mesh->u, mesh->dt);
    FILE *fp_u_star = fopen("../data/test_u_star.txt", "w+");
    save_header(t_u_star, mesh->u->x, mesh->u->y, fp_u_star);

    Tracker *t_v_star = track_mesh(mesh->v, mesh->dt);
    FILE *fp_v_star = fopen("../data/test_v_star.txt", "w+");
    save_header(t_v_star, mesh->v->x, mesh->v->y, fp_v_star);


    iterate(mesh, poisson, ic);

    // Save and close u*, v*
    save_tracking_state(t_u_star, 0.0, fp_u_star);
    fclose(fp_u_star);
    free_Tracker(t_u_star);

    save_tracking_state(t_v_star, 0.0, fp_v_star);
    fclose(fp_v_star);
    free_Tracker(t_v_star);

    printf("\t3.1) u*, v* trackers are saved and closed!\n");

    // Now that u* and v* are computed, we can compute the rhs
    printf("\t3.2) Tracking rhs\n");
    Tracker *t_rhs = track_mesh(mesh->p, mesh->dt);
    FILE *fp_rhs = fopen("../data/test_rhs.txt", "w+");
    save_header(t_rhs, mesh->p->x, mesh->p->y, fp_rhs);

    compute_rhs(mesh, mesh->p->val2, mesh->dt);

    save_tracking_state(t_rhs, 0.0, fp_rhs);
    fclose(fp_rhs);
    free_Tracker(t_rhs);

    printf("\t3.2) rhs tracker is saved and closed!\n");
}

int main(int argc, char *argv[]){


    PetscInitialize(&argc, &argv, 0, 0);

    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(CYLINDER);
    IterateCache *ic = initIterateCache(mesh);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh);

    // poisson_solver(poisson, mesh);
    // save_mesh(mesh->p);
    // printf("%d\n", mesh->p->n1);
    // printf("%d\n", mesh->p->n2);
    // exit(1);

    // for (int j = 0; j < mesh->u->n2; j++) {
    //     int ind_inner = j;
    //     int ind_outer = (mesh->u->n1-1)*(mesh->u->n2) + j;
    //
    //     mesh->u->val1[ind_outer] = cos(mesh->u->theta[ind_outer]);
    //     mesh->u->val1[ind_inner] = 0;
    // }

    int debug = 0;

    if (debug) {
        printf("Entered debugging mode\n");
        run_tests();
        exit(1);
    }


    double state = 0.0;         // time
    double dt    = mesh->dt;    // detla-time
    double nu    = NU;

    double endState = 1.0;

    int every_n = 200;
    int max_n = ceil(endState / dt);

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
                save_mesh_state(meshes[i], state * U_INF * 50.0, files[i]);
            }

        }

        /*
        if (ic->n == max_n) {
            break;
        }*/
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
