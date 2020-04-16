#include "solver.h"


/*
 *  Computes 1/dt * div(u*) at the mesh points (i,j)
 *  It is the right-hand side of the Poisson Equation.
 */
void compute_rhs(MACMesh *mesh, double *result, double dt) {
    double *u_star = mesh->u->val2;
    double *v_star = mesh->v->val2;

    int ind;
    double d1 = mesh->p->d1;
    double d2 = mesh->p->d2;
    double u_ij, v_ij;
    double du_d1, dv_d2;
    double h1, h2, dh1_d2, dh2_d1;
    for (int i = 0; i < mesh->p->n1; i++) {
        for (int j = 0; j < mesh->p->n2; j++) {
            ind = i*mesh->p->n2 + j;
            h1 = mesh->p->h1[ind];
            h2 = mesh->p->h2[ind];
            dh1_d2 = mesh->p->dh1_d2[ind];
            dh2_d1 = mesh->p->dh2_d1[ind];

            // Compute du*/dxi1
            du_d1 = (u_star[(i+1)*mesh->u->n2+j] - u_star[i*mesh->u->n2+j]) / (2*d1);

            // Compute dv*/dxi2
            if (j == mesh->v->n2-1) {   // We gotta pay attention to the periodicity
                dv_d2 = (v_star[i*mesh->v->n2]       - v_star[i*mesh->v->n2+j]) / (2*d2);
            } else {
                dv_d2 = (v_star[i*mesh->v->n2+(j+1)] - v_star[i*mesh->v->n2+j]) / (2*d2);
            }

            // Compute u_ij & v_ij
            u_ij = (u_star[i*mesh->u->n2+j] + u_star[(i+1)*mesh->u->n2+j] ) / 2.0;
            if (j == mesh->v->n2-1) {
                v_ij = (v_star[i*mesh->v->n2+j] + v_star[i*mesh->v->n2      ]) / 2.0;
            } else {
                v_ij = (v_star[i*mesh->v->n2+j] + v_star[i*mesh->v->n2+(j+1)]) / 2.0;
            }

            // Compute the complete term
            result[ind] = ((h2*du_d1 + dh2_d1*u_ij) + (h1*dv_d2 + dh1_d2*v_ij)) / (h1*h2*dt);
        }
    }

}


/*
 *  Computes the pressure gradient.
 *  On the border, we do not calculate those values: the velocities are known,
 *  therefore we do not need it.
 *  p_or_phi > 0 -> grad(pressure)
 *  p_or_phi < 0 -> grad(phi)
 */
void compute_grad_scalar(MACMesh *mesh, double *res_x, double *res_y, int p_or_phi) {
    double *field;
    if (p_or_phi > 0) {
        field = mesh->p->val1;
    } else {
        field = mesh->p->val2;
    }

    // Init some vars
    int ind;
    double d1 = mesh->p->d1;
    double d2 = mesh->p->d2;
    doublle h1, h2;

    // First, compute in the x-direction
    for (int i = 1; i < mesh->u->n1-1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            ind = i*mesh->u->n2 + j;
            h1 = mesh->u->h1[ind];
            h2 = mesh->u->h2[ind];

            res_x[ind] = (field[(i+1)*mesh->p->n2+j] - field[i*mesh->p->n2+j]) / (d1*h1);
        }
    }

    // Then compute in the y-direction
    for (int i = 1; i < mesh->v->n1-1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;
            h1 = mesh->v->h1[ind];
            h2 = mesh->v->h2[ind];

            if (j == mesh->v->n2-1) {
                res_y[ind] = (field[i*mesh->p->n2]       - field[i*mesh->p->n2+j]) / (d2*h2);
            } else {
                res_y[ind] = (field[i*mesh->p->n2+(j+1)] - field[i*mesh->p->n2+j]) / (d2*h2);
            }
        }
    }


}
