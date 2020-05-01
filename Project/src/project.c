#include "poisson.h"
#include "airfoil.h"
#include "cylinder.h"
#include "mesh.h"
#include "stdio.h"
#include "solver.h"
#include "debug.h"





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


    iterate(mesh, poisson, ic, 0);

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

    int debug = 0;

    if (debug) {
        printf("Entered debugging mode\n");
        run_tests();
        exit(1);
    }


    double state = 0.0;         // time
    double dt    = mesh->dt;    // detla-time
    double nu    = mesh->nu;

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
        iterate(mesh, poisson, ic, state);

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
