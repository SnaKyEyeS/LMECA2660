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

            // Between u(i+1/2, j) and u(i-1/2, j) there is d1 in xi1 distance,
            // so du(i, j) = 1/2 * (u(i+1/2, j) - u(i-1/2, j)) / (d1 / 2).
            // Same applies in xi2 direction.

            // Compute du*/dxi1
            du_d1 = (u_star[(i+1)*mesh->u->n2+j] - u_star[i*mesh->u->n2+j]) / d1;

            // Compute dv*/dxi2
            if (j == mesh->p->n2-1) {   // We gotta pay attention to the periodicity
                dv_d2 = (v_star[i*mesh->v->n2]       - v_star[i*mesh->v->n2+j]) / d2;
            } else {
                dv_d2 = (v_star[i*mesh->v->n2+(j+1)] - v_star[i*mesh->v->n2+j]) / d2;
            }

            // Compute u_ij & v_ij, the means of speeds are +-1/2 the current index
            u_ij = (u_star[i*mesh->u->n2+j] + u_star[(i+1)*mesh->u->n2+j] ) / 2.0;
            if (j == mesh->p->n2-1) {
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
    double h1, h2;

    // First, compute in the x-direction
    // We do not need the values at the boundaries Re and Ri, since those are
    // boundary conditions.
    for (int i = 1; i < mesh->u->n1-1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            ind = i*mesh->u->n2 + j;
            h1 = mesh->u->h1[ind];
            h2 = mesh->u->h2[ind];

            res_x[ind] = (field[(i+1)*mesh->p->n2+j] - field[i*mesh->p->n2+j]) / (d1*h1);
        }
    }

    // Then compute in the y-direction
    for (int i = 0; i < mesh->v->n1; i++) {
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


/*
 *  Computes w = rot(u).
 *  Evalutes it at points (i+1/2,j+1/2) and stores it in the associated mesh points.
 */
void compute_omega(MACMesh *mesh) {
    int ind;
    double d1 = mesh->w->d1;
    double d2 = mesh->w->d2;
    double h1, h2;
    double dh2_d1, dh1_d2;
    double dv_d1, du_d2;
    double *u = mesh->u->val1;
    double *v = mesh->v->val2;
    double v_ghost;

    for (int i = 0; i < mesh->w->n1; i++) {
        for (int j = 0; j < mesh->w->n2; j++) {
            ind = i*mesh->v->n2 + j;
            h1 = mesh->w->h1[ind];
            h2 = mesh->w->h2[ind];
            dh2_d1 = mesh->w->dh2_d1[ind];
            dh1_d2 = mesh->w->dh1_d2[ind];

            // 4th order finite differences
            // We use decentered schemes for the wall points.
            if (i == 0) {
                v_ghost = 0;    // COMPUTE THE VALUE using a polynomial
                dv_d1 = (-23*v_ghost                + 21*v[i*mesh->v->n2+j]    + 3 *v[(i+1)*mesh->v->n2+j] - 1*v[(i+2)*mesh->v->n2+j]) / (24*d1);
            } else if (i == 1) {
                dv_d1 = (-23*v[(i-1)*mesh->v->n2+j] + 21*v[i*mesh->v->n2+j]    + 3 *v[(i+1)*mesh->v->n2+j] - 1*v[(i+2)*mesh->v->n2+j]) / (24*d1);
            } else if (i == mesh->w->n1-2) {
                dv_d1 = (1*  v[(i-3)*mesh->v->n2+j] - 3*v[(i-2)*mesh->v->n2+j] - 21*v[(i-1)*mesh->v->n2+j] + 23*v[i*mesh->v->n2+j]) / (24*d1);
            } else if (i == mesh->w->n1-1) {
                v_ghost = 0;    // COMPUTE THE VALUE using a polynomial
                dv_d1 = (1*  v[(i-3)*mesh->v->n2+j] - 3*v[(i-2)*mesh->v->n2+j] - 21*v[(i-1)*mesh->v->n2+j] + 23*v_ghost           ) / (24*d1);
            } else {
                dv_d1 = (1*  v[(i-2)*mesh->v->n2+j] - 27*v[(i-1)*mesh->v->n2+j] + 27*v[i*mesh->v->n2+j] - 1*v[(i+1)*mesh->v->n2+j]) / (24*d1);
            }
            // We can use the periodicity in the xi2 direction (using the modulo operator)
            du_d2 = ( 1*u[i*mesh->u->n2+((j-2+mesh->u->n2)%mesh->u->n2)] - 27*u[i*mesh->u->n2+((j-1+mesh->u->n2)%mesh->u->n2)]
                   + 27*u[i*mesh->u->n2+(j%mesh->u->n2)]                  -  1*u[i*mesh->u->n2+((j+1)%mesh->u->n2)])              / (24*d2);

            mesh->w->val1[ind] = ((h2*dv_d1 + dh2_d1*v[ind]) - (h1*du_d2 + dh1_d2*u[ind])) / (h1*h2);
        }
    }
}


/*
 *  Computes the diffusive term nu*lapl(u).
 */
void compute_diffusive(MACMesh *mesh, double *res_x, double *res_y, double nu) {
    // We need to compute omega first
    compute_omega(mesh);

    // Init some vars
    int ind;
    double d1 = mesh->p->d1;
    double d2 = mesh->p->d2;
    double h1, h2;
    double *w = mesh->w->val1;

    // First, compute the x-composant --> the values at R_e and R_i are not needed
    // (boundary conditions).
    for (int i = 1; i < mesh->u->n1-1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            ind = i*mesh->u->n2 + j;
            h1 = mesh->u->h1[ind];
            h2 = mesh->u->h2[ind];

            if (j == mesh->w->n2) {
                res_x[ind] = -nu * (w[i*mesh->w->n2      ] - w[i*mesh->w->n2+j]) / (d2*h2);
            } else {
                res_x[ind] = -nu * (w[i*mesh->w->n2+(j+1)] - w[i*mesh->w->n2+j]) / (d2*h2);
            }
        }
    }

    // Then compute in the y-direction
    for (int i = 0; i < mesh->v->n1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;
            h1 = mesh->v->h1[ind];
            h2 = mesh->v->h2[ind];

            res_y[ind] = nu * (w[(i+1)*mesh->w->n2+j] - w[i*mesh->w->n2+j]) / (d1*h1);
        }
    }
}

// Compute h = u . grad u
// TODO !
void compute_h(MACMesh *mesh, double *res_x, double *res_y) {
    // Init some vars
    int ind;
    double d1 = mesh->d1;
    double d2 = mesh->d2;
    double h1, h2;
    double dh2_d1, dh1_d2;
    double *u = mesh->u->val1;
    double *v = mesh->v->val1;

    // First, compute the x-composant
    // Again, we do not need the values et Ri and Re since those are
    // boundary conditions.
    double du_d1, du_d2, v_avg;
    for (int i = 1; i < mesh->u->n1-1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            ind = i*mesh->u->n2 + j;
            h1 = mesh->u->h1[ind];
            h2 = mesh->u->h2[ind];
            dh2_d1 = mesh->u->dh2_d1[ind];
            dh1_d2 = mesh->u->dh1_d2[ind];

            // TODO/ periodic BC
            du_d1 = (u[ind+mesh->u->n2] - u[ind-mesh->u->n2]) / (2*d1);
            du_d2 = (u[ind+1]           - u[ind-1])           / (2*d2);
            v_avg = (v[(i-1)*mesh->v->n2+j] + v[(i-1)*mesh->v->n2+(j+1)]
                   + v[i*mesh->v->n2+(j+1)] + v[i*mesh->v->n2+j]         ) / 4;

            res_x[ind] = u[ind]*du_d1/h1 + v_avg*du_d2/h2 + v_avg*(u[ind]*dh1_d2 - v_avg*dh2_d1)/(h1*h2);
        }
    }

    // Then compute in the y-direction
    double dv_d1, dv_d2, u_avg;
    for (int i = 0; i < mesh->v->n1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;
            h1 = mesh->v->h1[ind];
            h2 = mesh->v->h2[ind];
            dh2_d1 = mesh->v->dh2_d1[ind];
            dh1_d2 = mesh->v->dh1_d2[ind];

            // TODO/ periodic BC
            dv_d1 = (v[ind+mesh->v->n2] - v[ind+mesh->v->n2]) / (2*d1);
            dv_d2 = (v[ind+1]           - v[ind-1])           / (2*d2);
            u_avg = (u[i*mesh->u->n2+(j-1)] + u[i*mesh->u->n2+j]
                    +u[(i+1)*mesh->u->n2+j] + u[(i+1)*mesh->u->n2+(j-1)]) / 4;

            res_y[ind] = u_avg*dv_d1/h1 + v[ind]*dv_d2/h2 + u_avg*(v[ind]*dh2_d1 - u_avg*dh1_d2)/(h1*h2);
        }
    }
}
