#include "solver.h"

// Compute indexes to handle periodicity
int index(int i, int j, int N, int i_shift, int j_shift) {
    return (i+i_shift) * N + ((N+j+j_shift) % N);
}


/*
 *  Computes 1/dt * div(u*) at the mesh points (i,j)
 *  It is the right-hand side of the Poisson Equation.
 */
void compute_rhs(MACMesh *mesh, double *result, double dt) {
    double *u_star = mesh->u->val2;
    double *v_star = mesh->v->val2;

    int ind;

    int ind_u_left, ind_u_right;
    int ind_v_bottom, ind_v_up;

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

            ind_u_left = index(i, j, mesh->u->n2, 0, 0);
            ind_u_right = index(i, j, mesh->u->n2, 1, 0);

            ind_v_bottom = index(i, j, mesh->v->n2, 0, 0);
            ind_v_up = index(i, j, mesh->v->n2, 0, 1);

            // Between u(i+1/2, j) and u(i-1/2, j) there is d1 in xi1 distance,
            // so du(i, j) = 1/2 * (u(i+1/2, j) - u(i-1/2, j)) / (d1 / 2).
            // Same applies in xi2 direction.

            // Compute du*/dxi1
            du_d1 = (u_star[ind_u_right] - u_star[ind_u_left]) / d1;

            // Compute dv*/dxi2
            dv_d2 = (v_star[ind_v_up] - v_star[ind_v_bottom]) / d2;

            // Compute u_ij & v_ij, the means of speeds are +-1/2 the current index
            u_ij = (u_star[ind_u_right] + u_star[ind_u_left])   / 2.0;
            v_ij = (v_star[ind_v_up]    + v_star[ind_v_bottom]) / 2.0;

            // Compute the complete term
            result[ind] = ((h2*du_d1 + dh2_d1*u_ij) + (h1*dv_d2 + dh1_d2*v_ij)) / (h1*h2*dt);
        }
    }
}


/*
 *  Computes the pressure gradient.
 *  On the border, we do not calculate those values: the velocities are known,
 *  therefore we do not need it.
 */
void compute_grad_scalar(MACMesh *mesh, double *res_x, double *res_y, GradientType type) {
    double *field;
    switch (type) {
        case PRESSURE:
            field = mesh->p->val1;
            break;

        case PHI:
            field = mesh->p->val2;
            break;

        default:
            printf("Gradient type not recognized !\n");
            exit(0);
    }

    // Init some vars
    int ind;
    int ind_p_left, ind_p_right, ind_p_bottom, ind_p_up;

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

            ind_p_left = index(i, j, mesh->p->n2, 0, 0);
            ind_p_right = index(i, j, mesh->p->n2, 1, 0);

            res_x[ind] = (field[ind_p_right] - field[ind_p_left]) / (d1*h1);
        }
    }

    // Then compute in the y-direction
    for (int i = 0; i < mesh->v->n1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;
            h1 = mesh->v->h1[ind];
            h2 = mesh->v->h2[ind];

            ind_p_bottom = index(i, j, mesh->p->n2, 0, 0);
            ind_p_up = index(i, j, mesh->p->n2, 0, 1);

            res_y[ind] = (field[ind_p_up] - field[ind_p_bottom]) / (d2*h2);
        }
    }
}


/*
 *  Computes w = rot(u).
 *  Evalutes it at points (i+1/2,j+1/2) and stores it in the associated mesh points.
 */
void compute_omega(MACMesh *mesh) {
    int ind;
    int ind_u_bottomx1, ind_u_bottomx2;
    int ind_u_upx1, ind_u_upx2;
    int ind_v_rightx1, ind_v_rightx2,  ind_v_rightx3, ind_v_rightx4;
    int ind_v_leftx1, ind_v_leftx2,  ind_v_leftx3, ind_v_leftx4;

    double d1 = mesh->w->d1;
    double d2 = mesh->w->d2;
    double h1, h2;
    double dh2_d1, dh1_d2;
    double dv_d1, du_d2;
    double *u = mesh->u->val1;
    double *v = mesh->v->val1;
    double v_ghost, v_avg, u_avg, v_wall, theta;

    for (int i = 0; i < mesh->w->n1; i++) {
        for (int j = 0; j < mesh->w->n2; j++) {
            ind = i*mesh->w->n2 + j;
            h1 = mesh->w->h1[ind];
            h2 = mesh->w->h2[ind];
            dh2_d1 = mesh->w->dh2_d1[ind];
            dh1_d2 = mesh->w->dh1_d2[ind];

            ind_u_upx2 = index(i, j, mesh->u->n2, 0, 1);
            ind_u_upx1 = index(i, j, mesh->u->n2, 0, 0);
            ind_u_bottomx1 = index(i, j, mesh->u->n2, 0, -1);
            ind_u_bottomx2 = index(i, j, mesh->u->n2, 0, -2);

            ind_v_rightx3 = index(i, j, mesh->v->n2, 2, 0);
            ind_v_rightx2 = index(i, j, mesh->v->n2, 1, 0);
            ind_v_rightx1 = index(i, j, mesh->v->n2, 0, 0);
            ind_v_leftx1 = index(i, j, mesh->v->n2, -1, 0);
            ind_v_leftx2 = index(i, j, mesh->v->n2, -2, 0);
            ind_v_leftx3 = index(i, j, mesh->v->n2, -3, 0);

            // 4th order finite differences
            // We use decentered schemes for the wall points
            // and a ghost point whose value is dictated by the BC for v.
            if (i == 0) {
                // v_wall = 0;
                // v_ghost = (5*v[ind_v_rightx4] - 28*v[ind_v_rightx3] + 70*v[ind_v_rightx2] - 140*v[ind_v_rightx1] + 128*v_wall)/35;
                // dv_d1 = (-23*v_ghost         + 21*v[ind_v_rightx1]    + 3 *v[ind_v_rightx2] - 1*v[ind_v_rightx3]) / (24*d1);
                ind_v_rightx4 = index(i, j, mesh->v->n2, 3, 0);
                dv_d1 = (-71*v[ind_v_rightx1] + 141*v[ind_v_rightx2] - 93*ind[ind_v_rightx3] + 23*v[ind_v_rightx4]) / (24*d1);
            } else if (i == 1) {
                dv_d1 = (-23*v[ind_v_leftx1] + 21*v[ind_v_rightx1]    + 3 *v[ind_v_rightx2] - 1*v[ind_v_rightx3]) / (24*d1);

            } else if (i == mesh->w->n1-2) {
                dv_d1 = (1*  v[ind_v_leftx3] - 3*v[ind_v_leftx2] - 21*v[ind_v_leftx1] + 23*v[ind_v_rightx1]) / (24*d1);
            } else if (i == mesh->w->n1-1) {
                // ind_v_leftx4 = index(i, j, mesh->v->n2, -4, 0);
                // theta = mesh->w->theta[ind];
                // v_wall = -U_INF * sin(theta) + U_PERT * cos(theta); // TODO: valeur de v_wall en r = R_ext ????
                // v_ghost = (5*v[ind_v_leftx4] - 28*v[ind_v_leftx3] + 70*v[ind_v_leftx2] - 140*v[ind_v_leftx1] + 128*v_wall)/35;
                // dv_d1 = (1*  v[ind_v_leftx3] - 3*v[ind_v_leftx2] - 21*v[ind_v_leftx1] + 23*v_ghost           ) / (24*d1);
                
            } else {
                dv_d1 = (1*  v[ind_v_leftx2] - 27*v[ind_v_leftx1] + 27*v[ind_v_rightx1] - 1*v[ind_v_rightx2]) / (24*d1);
            }
            // We can use the periodicity in the xi2 direction (using the modulo operator)
            // NB: the boundary conditions at the walls must be imposed in u[] before the call to this function !
            du_d2 = ( 1*u[ind_u_bottomx2] - 27*u[ind_u_bottomx1] + 27*u[ind_u_upx1] - 1*u[ind_u_upx2]) / (24*d2);

            // Compute v & u averaged
            u_avg = (u[ind_u_upx1]    + u[ind_u_bottomx1]) / 2;
            if (i == 0) {
                v_avg = 0;
            } else if (i == mesh->w->n1) {
                theta = mesh->w->theta[ind];
                v_avg = -U_INF * sin(theta) + U_PERT * cos(theta);
                            // TODO: trouver quelle est la valeur de v sur le mur, je sais pas ce que ça vaut :(
                            // doit être calculée telle que w_wall = 0...
            } else {
                v_avg = (v[ind_v_rightx1] + v[ind_v_leftx1])   / 2;
            }

            if (i == mesh->w->n1-1) {
                mesh->w->val1[ind] = 0.0;   // Set w = 0 at the outerwall
            } else {
                mesh->w->val1[ind] = ((h2*dv_d1 + dh2_d1*v_avg) - (h1*du_d2 + dh1_d2*u_avg)) / (h1*h2);
            }
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

    int ind_w_left, ind_w_right, ind_w_bottom, ind_w_up;

    double d1 = mesh->p->d1;
    double d2 = mesh->p->d2;
    double h1, h2;
    double *w = mesh->w->val1;

    // First, compute the x-composant --> the values at R_e and R_i are not needed
    // (boundary conditions).
    for (int i = 0; i < mesh->u->n1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            ind = i*mesh->u->n2 + j;
            h1 = mesh->u->h1[ind];
            h2 = mesh->u->h2[ind];

            ind_w_bottom = index(i, j, mesh->w->n2, 0, 0);
            ind_w_up = index(i, j, mesh->w->n2, 0, 1);

            res_x[ind] = -nu * (w[ind_w_up] - w[ind_w_bottom]) / (d2*h2);
        }
    }

    // Then compute in the y-direction
    for (int i = 0; i < mesh->v->n1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;
            h1 = mesh->v->h1[ind];
            h2 = mesh->v->h2[ind];

            ind_w_left = index(i, j, mesh->w->n2, 0, 0);
            ind_w_right = index(i, j, mesh->w->n2, 1, 0);

            res_y[ind] = nu * (w[ind_w_right] - w[ind_w_left]) / (d1*h1);
        }
    }
}

// Compute h = u . grad u
void compute_h(MACMesh *mesh, double *res_x, double *res_y) {
    // Init some vars
    int ind;
    int ind_u_left, ind_u_right, ind_u_up, ind_u_bottom;
    int ind_v_bottom_left, ind_v_bottom_right, ind_v_up_left, ind_v_up_right;
    double theta, u_wall;
    double d1 = mesh->d1;
    double d2 = mesh->d2;
    double h1, h2;
    double dh2_d1, dh1_d2;
    double *u = mesh->u->val1;
    double *v = mesh->v->val1;

    // First, compute the x-composant
    // Again, we do not need the values at Ri and Re since those are
    // boundary conditions, thus are imposed.
    double du_d1, du_d2, v_avg;
    for (int i = 1; i < mesh->u->n1-1; i++) {
        for (int j = 0; j < mesh->u->n2; j++) {
            // We compute the indexes
            ind = i*mesh->u->n2 + j;

            ind_u_left          = index(i, j, mesh->u->n2, -1, 0);
            ind_u_right         = index(i, j, mesh->u->n2, 1, 0);
            ind_u_up            = index(i, j, mesh->u->n2, 0, 1);
            ind_u_bottom        = index(i, j, mesh->u->n2, 0, -1);

            ind_v_bottom_left   = index(i, j, mesh->v->n2, -1, 0);
            ind_v_bottom_right  = index(i, j, mesh->v->n2, 0, 0);
            ind_v_up_left       = index(i, j, mesh->v->n2, -1, 1);
            ind_v_up_right      = index(i, j, mesh->v->n2, 0, 1);

            h1 = mesh->u->h1[ind];
            h2 = mesh->u->h2[ind];
            dh2_d1 = mesh->u->dh2_d1[ind];
            dh1_d2 = mesh->u->dh1_d2[ind];

            du_d1 = (u[ind_u_right] - u[ind_u_left])   / (2*d1);
            du_d2 = (u[ind_u_up]    - u[ind_u_bottom]) / (2*d2);
            v_avg = (v[ind_v_bottom_left]  + v[ind_v_bottom_right]
                    +v[ind_v_up_left]      + v[ind_v_up_right]    ) / 4;

            res_x[ind] = u[ind]*du_d1/h1 + v_avg*du_d2/h2 + v_avg*(u[ind]*dh1_d2 - v_avg*dh2_d1)/(h1*h2);
        }
    }

    int ind_w_wall;
    int ind_v_right_right, ind_v_left_left;
    int ind_v_left, ind_v_right, ind_v_up, ind_v_bottom;
    int ind_u_bottom_left, ind_u_bottom_right, ind_u_up_left, ind_u_up_right;
    double theta_u_up_right, theta_u_bottom_right;
    double u_wall_up_right, u_wall_bottom_right;
    double v_ghost_left, v_ghost_right;
    double v_wall_right;

    // Then compute in the y-direction
    double dv_d1, dv_d2, u_avg;
    for (int i = 0; i < mesh->v->n1; i++) {
        for (int j = 0; j < mesh->v->n2; j++) {
            ind = i*mesh->v->n2 + j;

            ind_v_left          = index(i, j, mesh->v->n2, -1, 0);
            ind_v_right         = index(i, j, mesh->v->n2, 1, 0);
            ind_v_up            = index(i, j, mesh->v->n2, 0, 1);
            ind_v_bottom        = index(i, j, mesh->v->n2, 0, -1);

            ind_u_bottom_left   = index(i, j, mesh->u->n2, 0, -1);
            ind_u_bottom_right  = index(i, j, mesh->u->n2, 1, -1);
            ind_u_up_left       = index(i, j, mesh->u->n2, 0, 0);
            ind_u_up_right      = index(i, j, mesh->u->n2, 1, 0);

            h1 = mesh->v->h1[ind];
            h2 = mesh->v->h2[ind];
            dh2_d1 = mesh->v->dh2_d1[ind];
            dh1_d2 = mesh->v->dh1_d2[ind];

            if (i == 0) {                                   // If at r = Ri
                ind_v_right_right = index(i, j, mesh->v->n2, 2, 0);
                v_ghost_left = -(v[ind_v_right_right] - 5*v[ind_v_right] + 15*v[ind] - 16*0)/5;   // We want v_wall = 0 !
                dv_d1 = (v[ind_v_right] - v_ghost_left)    / (2*d1);
                dv_d2 = (v[ind_v_up]    - v[ind_v_bottom]) / (2*d2);
                u_avg = (   0.0      + u[ind_u_bottom_right]
                        +   0.0      + u[ind_u_up_right]    ) / 4;
            }
            // need to check the interpolation (I took the opposite of that of the left wall)
            // and the wall velocity for v !
            else if (i == mesh->v->n1-1) {                      // If at r = Re
                ind_w_wall = index(i, j, mesh->w->n2, 1, 0);    // We need the orientation of the mesh at the wall
                ind_v_left_left = index(i, j, mesh->v->n2, -2, 0);
                theta = mesh->w->theta[ind_w_wall];
                v_wall_right = -U_INF * sin(theta) + U_PERT * cos(theta);       // TODO valeur de v_wall en r = Ri... j'suis pas sur de ce truc
                v_ghost_right = -(v[ind_v_left_left] - 5*v[ind_v_left] + 15*v[ind] - 16*v_wall_right)/5;
                dv_d1 = (v_ghost_right - v[ind_v_left])   / (2*d1);
                dv_d2 = (v[ind_v_up]   - v[ind_v_bottom]) / (2*d2);
                u_avg = (   u[ind_u_bottom_left]  + u[ind_u_bottom_right]
                        +   u[ind_u_up_left]      + u[ind_u_up_right]    ) / 4;
            }
            else {
                dv_d1 = (v[ind_v_right] - v[ind_v_left]) / (2*d1);
                dv_d2 = (v[ind_v_up]           - v[ind_v_bottom])           / (2*d2);
                u_avg = (   u[ind_u_bottom_left]  + u[ind_u_bottom_right]
                        +   u[ind_u_up_left]      + u[ind_u_up_right]    ) / 4;
            }

            res_y[ind] = u_avg*dv_d1/h1 + v[ind]*dv_d2/h2 + u_avg*(v[ind]*dh2_d1 - u_avg*dh1_d2)/(h1*h2);
        }
    }
}
