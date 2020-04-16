#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "mesh.h"


// Function declaration
void compute_rhs(MACMesh *mesh, double *result, double dt);
void compute_grad_scalar(MACMesh *mesh, double *res_x, double *res_y, int p_or_phi);


#endif
