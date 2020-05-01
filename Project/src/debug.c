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
