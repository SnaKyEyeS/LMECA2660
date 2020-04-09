#include "cylinder.h"


CylinderMapping *init_cylinder_mapping() {
    CylinderMapping *mapping = (CylinderMapping *) malloc(sizeof(CylinderMapping));

    // Fixed parameters
    mapping->D = 0.2;                           // Fixé arbitrairement ???
    mapping->R = mapping->D / 2.0;
    mapping->H = 50*mapping->D;                 // H = 50D, OK
    mapping->h_wall_normal = 2.0e-4;            // On doit déterminer ces param pour avoir
    mapping->gamma = 1.015;                     // Re_w = abs(w)*h_max² / nu < 40 dans la région 12D.
    mapping->n_xi2 = 360;

    // We compute the rest
    mapping->beta = mapping->h_wall_normal / (mapping->gamma - 1);
    mapping->alpha = log(mapping->H/mapping->beta + 1);
    mapping->dxi1 = log(mapping->gamma) / mapping->alpha;
    mapping->n_xi1 = ceil(1 / mapping->dxi1);
    mapping->dxi2 = 2*M_PI / mapping->n_xi2;        // Pareil ici ???

    // Mapping limits
    mapping->xi1_lim[0] = 0.0;
    mapping->xi1_lim[1] = 1.0;

    mapping->xi2_lim[0] = 0.0;
    mapping->xi2_lim[1] = 2*M_PI;

    printf("%d, %d\n", mapping->n_xi1);
    return mapping;
}


void cylinder_metrics(CylinderMapping *mapping, double xi1, double xi2, double *x, double *y,
                      double *h1, double *h2, double *dh1_dxi1, double *dh1_dxi2, double *dh2_dxi1, double *dh2_dxi2,
                      double *d2h1_dxi1dxi2, double *d2h2_dxi1dxi2, double *theta_mesh) {
    // X and Y coordinates
    if (x) { *x = (mapping->R + mapping->beta*(exp(xi1*mapping->alpha) - 1)) * cos(xi2); }
    if (y) { *y = (mapping->R + mapping->beta*(exp(xi1*mapping->alpha) - 1)) * sin(xi2); }

    // Form factor H1 & H2
    if (h1) { *h1 = SQRT_2 * mapping->alpha * mapping->beta * exp(xi1 * mapping->alpha); }
    if (h2) { *h2 = SQRT_2 * (mapping->R + mapping->beta * (exp(xi1 * mapping->alpha) - 1)); }

    // First derivative of the form factor
    if (dh1_dxi1) { *dh1_dxi1 = SQRT_2 * mapping->alpha*mapping->alpha * mapping->beta * exp(xi1 * mapping->alpha); }
    if (dh1_dxi2) { *dh1_dxi2 = 0.0; }
    if (dh2_dxi1) { *dh2_dxi1 = SQRT_2 * mapping->alpha*mapping->alpha * mapping->beta * exp(xi1 * mapping->alpha); }
    if (dh2_dxi2) { *dh2_dxi2 = 0.0; }

    // Second derivative of the form factor
    if (d2h2_dxi1dxi2) { *d2h2_dxi1dxi2 = 0.0; }
    if (d2h1_dxi1dxi2) { *d2h1_dxi1dxi2 = 0.0; }

    // Mesh orientation
    if (theta_mesh) { *theta_mesh = xi2 + HALF_PI; }
}

void cylinder_init_velocity(CylinderMapping *mapping, double U_inf, double xi1, double xi2, double *u_n, double *u_theta) {
    // ...
    double theta;
    double r = mapping->R + mapping->beta*(exp(xi1*mapping->alpha) - 1);
    cylinder_metrics(mapping, xi1, xi2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &theta);
    *u_n =       U_inf * (1 - (mapping->R*mapping->R) / (r*r)) * cos(theta);
    *u_theta = - U_inf * (1 + (mapping->R*mapping->R) / (r*r)) * sin(theta);
}

void cylinder_init_pressure(CylinderMapping *mapping, double U_inf, double xi1, double xi2, double *p) {
    // ...
    double u_n, u_t;
    cylinder_init_velocity(mapping, U_inf, xi1, xi2, &u_n, &u_t);
    *p = - fabs(u_n*u_n + u_t*u_t) / 2.0;
}
