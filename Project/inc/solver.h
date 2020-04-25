#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "mesh.h"

typedef enum {PRESSURE, PHI} GradientType;

// Function declaration
int index(int i, int j, int N, int i_shift, int j_shift);
void compute_rhs(MACMesh *mesh, double *result, double dt);
void compute_grad(MACMesh *mesh, double *res_x, double *res_y, GradientType type);
void compute_omega(MACMesh *mesh);
void compute_diffusive(MACMesh *mesh, double *res_x, double *res_y, double nu);
double interpolate2D(double x_1, double x_2, double y_1, double y_2, double U[4], double x, double y);
void compute_h(MACMesh *mesh, double *res_x, double *res_y);
#endif
