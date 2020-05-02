#include "cylinder.h"


CylinderMapping *init_cylinder_mapping() {
    CylinderMapping *mapping = (CylinderMapping *) malloc(sizeof(CylinderMapping));

    // Fixed parameters
    mapping->D = 1.0/50.0;
    mapping->R = mapping->D / 2.0;
    mapping->H = 50*mapping->D;
    mapping->h_wall_normal = mapping->D / 300;
    mapping->gamma = 1.015;
    mapping->n_xi2 = 360;

    // We compute the rest
    mapping->beta = mapping->h_wall_normal / (mapping->gamma - 1);
    mapping->alpha = log(mapping->H/mapping->beta + 1);
    mapping->dxi1 = log(mapping->gamma) / mapping->alpha;
    mapping->n_xi1 = ceil(1 / mapping->dxi1);
    mapping->dxi2 = 2*M_PI / mapping->n_xi2;        // Pareil ici ???

    // Mapping limits
    mapping->xi1_lim[0] = 0.0;
    mapping->xi1_lim[1] = mapping->n_xi1*mapping->dxi1;

    mapping->xi2_lim[0] = 0.0;
    mapping->xi2_lim[1] = 2*M_PI;


    mapping->Lc = mapping->D;
    return mapping;
}


void cylinder_metrics(CylinderMapping *mapping, double xi1, double xi2, double *x, double *y,
                      double *h1, double *h2, double *dh1_dxi1, double *dh1_dxi2, double *dh2_dxi1, double *dh2_dxi2,
                      double *d2h1_dxi1dxi2, double *d2h2_dxi1dxi2, double *theta_mesh) {

    double a = mapping->alpha;
    double b = mapping->beta;
    double R = mapping->R;

    // X and Y coordinates
    if (x) { *x = (R + b*(exp(xi1*a) - 1)) * cos(xi2); }
    if (y) { *y = (R + b*(exp(xi1*a) - 1)) * sin(xi2); }

    // Form factor H1 & H2
    if (h1) { *h1 = a * b * exp(xi1 * a); }
    if (h2) { *h2 = R + b * (exp(xi1 * a) - 1); }

    // First derivative of the form factor
    if (dh1_dxi1) { *dh1_dxi1 = a * a * b * exp(xi1 * a); }
    if (dh1_dxi2) { *dh1_dxi2 = 0.0; }
    if (dh2_dxi1) { *dh2_dxi1 = a * b * exp(xi1 * a); }
    if (dh2_dxi2) { *dh2_dxi2 = 0.0; }

    // Second derivative of the form factor
    if (d2h2_dxi1dxi2) { *d2h2_dxi1dxi2 = 0.0; }
    if (d2h1_dxi1dxi2) { *d2h1_dxi1dxi2 = 0.0; }

    // Mesh orientation
    if (theta_mesh) { *theta_mesh = xi2; }
}

void cylinder_init_velocity(CylinderMapping *mapping, double U_inf, double xi1, double xi2, double *u_n, double *u_t) {
    // ...
    double theta;
    double r = mapping->R + mapping->beta*(exp(xi1*mapping->alpha) - 1);
    cylinder_metrics(mapping, xi1, xi2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &theta);
    if(u_n) { *u_n =   U_inf * (1 - (mapping->R*mapping->R) / (r*r)) * cos(theta); }
    if(u_t) { *u_t = - U_inf * (1 + (mapping->R*mapping->R) / (r*r)) * sin(theta); }
}

void cylinder_init_pressure(CylinderMapping *mapping, double U_inf, double xi1, double xi2, double *p) {
    // ...
    double u_n, u_t;
    cylinder_init_velocity(mapping, U_inf, xi1, xi2, &u_n, &u_t);
    *p = - (u_n*u_n + u_t*u_t) / 2.0;
}
