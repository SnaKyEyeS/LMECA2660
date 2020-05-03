#include "debug.h"

/*
 * Pre-computed analytical values for cylinder statement
 * i.e. div u = 0
 * In this particular case of u = [(1 - R^2/r^2) * cos(theta),  - (1 + R^2/r^2) * sin(theta)] :
 *                    w = rot u = [0, 0, 3/2 * cos(theta)]
 *                  laplacian u = [- (cos(theta) * (r^2 + 3 * R^2)) / r^4, (sin(theta) * (r^2 - 3 * R^2)) / r^4]
 *                       grad u = [1/2 * sin(theta), -1/2 * cos(theta);
 *                                 cos(theta)      , -1/2 * sin(theta)]
 *             H = u dot grad u = [r/4 * sin^2(theta) - r/2 * cos^2(theta), sin(theta) * cos(theta) * (r/2 - r/2)]
 *                              = [r/4 * sin^2(theta) - r/2 * cos^2(theta), 0]
 *                            P = - |u|^2 / 2= - 1/8 * r^2 * sin^2(theta)
 *                       grad P = [- 1/4 * r * sin^2(theta), - 1/4 * r * sin(theta) * cos(theta)] // Higher error on y (?)
 * (first iteration)         u* = u + dt * (-H - grad P + nu * laplacian u)
 *                          rhs = div . u* / dt
 *   (div . u) / dt:            = 0 / dt
 * - (div . H):                   + 1/4 * (cos(2 * theta) + 3)
 * - (div . grad P):              + 1/4
 *   (div . laplacian u):         + 0
 *                              = 1/4 * cos(2 * theta) + 1
 *
 */

/*
 * Fills the speed field in cylindrical coordinates with a solenoidal function
 * i.e. div u = 0
 * In this particular case of u = [1/2 * r * sin(theta), r * cos(theta)] :
 *                    w = rot u = [0, 0, 3/2 * cos(theta)]
 *                  laplacian u = [3/2 * sin(theta) / r, 0, 0]
 *                       grad u = [1/2 * sin(theta), -1/2 * cos(theta);
 *                                 cos(theta)      , -1/2 * sin(theta)]
 *             H = u dot grad u = [r/4 * sin^2(theta) - r/2 * cos^2(theta), sin(theta) * cos(theta) * (r/2 - r/2)]
 *                              = [r/4 * sin^2(theta) - r/2 * cos^2(theta), 0]
 *                            P = - |u|^2 / 2= - 1/8 * r^2 * sin^2(theta)
 *                       grad P = [- 1/4 * r * sin^2(theta), - 1/4 * r * sin(theta) * cos(theta)] // Higher error on y (?)
 * (first iteration)         u* = u + dt * (-H - grad P + nu * laplacian u)
 *                          rhs = div . u* / dt
 *   (div . u) / dt:            = 0 / dt
 * - (div . H):                   + 1/4 * (cos(2 * theta) + 3)
 * - (div . grad P):              + 1/4
 *   (div . laplacian u):         + 0
 *                              = 1/4 * cos(2 * theta) + 1
 *
 */
void solenoidal_speed(double r, double theta, double *u, double *v) {
    if (u) {*u = 0.5 * r * sin(theta);}
    if (v) {*v =       r * cos(theta);}
}

/*
 * Tracker will track values in both normal and tangential directions
 */
Tracker *init_Tracker(double n1, double n2, double *val1, double *val2, double dt) {
    Tracker *t = (Tracker *) malloc(sizeof(Tracker));
    t->n1 = n1;
    t->n2 = n2;
    t->val1 = val1;
    t->val2 = val2;
    t->n = n1 * n2;
    t->dt = dt;

    return t;
}

void free_Tracker(Tracker *t) {
    free(t);
}

Tracker *track_mesh(Mesh *mesh, double dt) {
    return init_Tracker(mesh->n1, mesh->n2, mesh->val1, mesh->val2, dt);
}

void save_header(Tracker *t, double *x, double *y, FILE *fp) {
    fprintf(fp, "{'n1': %d, 'n2': %d, 'dt':%.10f, 'nu':%.10f}\n", t->n1, t->n2, t->dt, NU);
    save_array(x, t->n, fp);
    save_array(y, t->n, fp);
}

void save_tracking_state(Tracker *t, double state, FILE *fp) {
    fprintf(fp, "%.10f, ", state);
    save_array(t->val1, t->n, fp);
    fprintf(fp, "%.10f, ", state);
    save_array(t->val2, t->n, fp);
}

void set_speed_pressure(MACMesh *mesh, void (*f)(double, double, double*, double*)) {
    double x, y, r, theta;
    for (int ind = 0; ind < mesh->u->n; ind++) {
        x = mesh->u->x[ind];
        y = mesh->u->y[ind];
        r = hypot(x, y);
        theta = mesh->u->theta[ind];

        (*f)(r, theta, &mesh->u->val1[ind], NULL);
    }
    for (int ind = 0; ind < mesh->v->n; ind++) {
        x = mesh->v->x[ind];
        y = mesh->v->y[ind];
        r = hypot(x, y);
        theta = mesh->v->theta[ind];

        (*f)(r, theta, NULL, &mesh->v->val1[ind]);
    }
    // Presure from speed
    double u;
    for (int ind = 0; ind < mesh->p->n; ind++) {
        x = mesh->p->x[ind];
        y = mesh->p->y[ind];
        r = hypot(x, y);
        theta = mesh->p->theta[ind];

        (*f)(r, theta, &u, NULL);

        mesh->p->val1[ind] = -(u * u) * 0.5;
    }
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
