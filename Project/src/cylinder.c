#include "cylinder.h"


CylinderMapping *init_cylinder_mapping() {
    CylinderMapping *mapping = (CylinderMapping *) malloc(sizeof(CylinderMapping));

    // Fixed parameters
    mapping->D = 2;                           // Fixé arbitrairement ???
    mapping->R = mapping->D / 2.0;
    mapping->H = 20; //50*mapping->D;                 // H = 50D, OK
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

    printf("%d\n", mapping->n_xi1);
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
    if (theta_mesh) {
        double dxdxi1 = a * b * exp(xi1*a) * cos(xi2);
        double dydxi1 = a * b * exp(xi1*a) * sin(xi2);
        *theta_mesh = atan2(dydxi1, dxdxi1);        // Verifier sa definition ?
    }
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
    *p = - fabs(u_n*u_n + u_t*u_t) / 2.0;
}
