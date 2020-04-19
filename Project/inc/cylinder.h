#ifndef _CYLINDER_H_
#define _CYLINDER_H_

#include <math.h>
#include "const.h"
#include "stdio.h"

typedef struct {
    // Number of points
    int n_xi1;
    int n_xi2;

    // Discretization
    double dxi1;
    double dxi2;

    // Parameters of the Mapping
    double alpha;
    double beta;
    double gamma;
    double h_wall_normal;

    // Diameter of the cylinder
    double D;
    double R;

    // Distance between the cylinder and the outer boundary
    double H;

    // Mapping limits
    double xi1_lim[2];
    double xi2_lim[2];

    // Caracteristic length (for Reynolds numbers and such)
    double Lc;

} CylinderMapping;


// Function declarations
CylinderMapping*    init_cylinder_mapping();
void cylinder_metrics(CylinderMapping *mapping, double xi1, double xi2, double *x, double *y,
                                      double *h1, double *h2, double *dh1_dxi1, double *dh1_dxi2, double *dh2_dxi1, double *dh2_dxi2,
                                      double *d2h1_dxi1dxi2, double *d2h2_dxi1dxi2, double *theta_mesh);
void cylinder_init_velocity(CylinderMapping *mapping, double U_inf, double xi1, double xi2, double *u_n, double *u_theta);
void cylinder_init_pressure(CylinderMapping *mapping, double U_inf, double xi1, double xi2, double *p);

#endif
