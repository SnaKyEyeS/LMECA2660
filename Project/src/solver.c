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
