#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "thomas.h"

/*
    Choice of the scheme for the convective term. One can choose from:
     - E2 (2nd order explicit)
     - E4 (4th order explicit)
     - DS (3rd order decentered explicit)
     - I4 (4th order implicit)
     - I6 (6th order implicit)
*/
typedef enum {E2, E4, DS, I4, I6} Scheme;
#define FD_SCHEME I6

// Physical constants & mesh discretization
#define LENGTH          1.0
#define T_END           1.0
#define SIGMA_0         (LENGTH/32.0)
#define CONV_COEFF      (LENGTH/T_END)
#define DIFF_COEFF      (CONV_COEFF*SIGMA_0/40.0)        // 0 or (CONV_COEFF*SIGMA_0/40.0)
#define TIME_DISCR      (SPACE_DISCR/CONV_COEFF)
#define SPACE_DISCR     (SIGMA_0*0.125)                   // 0.5, 0.25 or 0.125
#define N_DISCR         (int)(LENGTH/SPACE_DISCR)


// Struct definitions
typedef struct Mesh {
    int n;
    double h;
    double *U;
} Mesh;

// Numerical schemes
void diffusive_term(double *U, int n, double h, double *d2U);
void convective_term(double *U, int n, double h, double *dU);
void runge_kutta_4(Mesh *mesh, double dt);
void compute_global_diagnostics(Mesh *mesh, double t, double *u, double *Q, double *E, double *R);
void compute_solution(double t, int n, int n_step, double *u, double *e);

// Init/free routines
Mesh *init_mesh();
void free_mesh(Mesh *mesh);

#endif
