#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "mesh.h"
#include "poisson.h"

typedef enum {PRESSURE, PHI} GradientType;

typedef struct {
    int n;          // Current iteration

    // Cache
    double *new_h_x, *new_h_y;
    double *old_h_x, *old_h_y;
    double *grad_P_x, *grad_P_y;
    double *grad_Phi_x, *grad_Phi_y;
    double *nu_lapl_u_x, *nu_lapl_u_y;
} IterateCache;


// Scheme discretisations
int index(int i, int j, int N, int i_shift, int j_shift);
void compute_rhs(MACMesh *mesh, double *result, double dt);
void compute_grad(MACMesh *mesh, double *res_x, double *res_y, GradientType type);
void compute_omega(MACMesh *mesh);
void compute_diffusive(MACMesh *mesh, double *res_x, double *res_y, double nu);
double interpolate2D(double x_1, double x_2, double y_1, double y_2, double U[4], double x, double y);
void compute_h(MACMesh *mesh, double *res_x, double *res_y);

// Iteration Cache
IterateCache *initIterateCache(MACMesh *mesh);
void freeIterateCache(IterateCache *ic);
void iterate(MACMesh *mesh, PoissonData *poisson, IterateCache *ic);
#endif
