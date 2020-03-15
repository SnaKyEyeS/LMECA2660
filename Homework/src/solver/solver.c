#include "solver.h"


void diffusive_term(double *U, int n, double h, double *d2U) {
    // Compute the second derivative (E2 implicit scheme)
    for (int i = 0; i < n; i++) {
        d2U[i] = (U[(i+1) % n] - 2*U[i] + U[(n+i-1) % n]) / (h*h);
    }
}

void convective_term(double *U, int n, double h, double *dU) {
    // Compute the first derivative
    switch (FD_SCHEME) {
        case E2:
            for (int i = 0; i < n; i++) {
                dU[i] = (U[(i+1) % n] - U[(n+i-1) % n]) / (2*h);
            }
            break;

        case E4:
            for (int i = 0; i < n; i++) {
                dU[i] = 4.0/3.0 * (U[(i+1) % n] - U[(n+i-1) % n]) / (2*h) - 1.0/3.0 * (U[(i+2) % n] - U[(n+i-2) % n]) / (4*h);
            }
            break;

        case DS:
            for (int i = 0; i < n; i++) {
                dU[i] = (U[(i+1) % n]/3.0 - U[(n+i-1) % n] + U[(n+i-2) % n]/6.0 + U[i]/2.0) / h;
            }
            break;

        case I4: {
            double *q = (double*) calloc(n, sizeof(double));
            for (int i = 0; i < n; i++) {
                q[i] = 3.0/2.0 * (U[(i+1) % n] - U[(n+i-1) % n]) / (2*h);
            }
            solve_period_3diag(n, 1.0, 0.25, 0.25, dU, q);
            free(q);
            break;
        }

        case I6: {
            double *q = (double*) calloc(n, sizeof(double));
            for (int i = 0; i < n; i++) {
                q[i] = 14.0/9.0 * (U[(i+1) % n] - U[(n+i-1) % n]) / (2*h) + 1.0/9.0 * (U[(i+2) % n] - U[(n+i-2) % n]) / (4*h);
            }
            solve_period_3diag(n, 1.0, 1.0/3.0, 1.0/3.0, dU, q);
            free(q);
            break;
        }

        default:
            printf("The requested scheme is not implemented. Please specify another one.");
            exit(0);
    }
}

void runge_kutta_4(Mesh *mesh, double dt) {
    // Initialize stuff
    int i;
    int n = mesh->n;
    double h = mesh->h;
    double *dU  = (double*) calloc(n, sizeof(double));
    double *d2U = (double*) calloc(n, sizeof(double));
    double *Up  = (double*) calloc(n, sizeof(double));
    double *K1  = (double*) calloc(n, sizeof(double));
    double *K2  = (double*) calloc(n, sizeof(double));
    double *K3  = (double*) calloc(n, sizeof(double));
    double *K4  = (double*) calloc(n, sizeof(double));

    // Compute K1, K2, K3 and K4
    // K1:
    diffusive_term(mesh->U, n, h, d2U);
    convective_term(mesh->U, n, h, dU);
    for(i = 0; i < n; i++) {
        K1[i] = DIFF_COEFF*d2U[i] - CONV_COEFF*dU[i];
    }

    // K2:
    for(i = 0; i < n; i++) {
        Up[i] = mesh->U[i] + K1[i] * dt/2.0;
    }
    diffusive_term(Up, n, h, d2U);
    convective_term(Up, n, h, dU);
    for(i = 0; i < n; i++) {
        K2[i] = DIFF_COEFF*d2U[i] - CONV_COEFF*dU[i];
    }

    // K3:
    for(i = 0; i < n; i++) {
        Up[i] = mesh->U[i] + K2[i] * dt/2.0;
    }
    diffusive_term(Up, n, h, d2U);
    convective_term(Up, n, h, dU);
    for(i = 0; i < n; i++) {
        K3[i] = DIFF_COEFF*d2U[i] - CONV_COEFF*dU[i];
    }

    // K4:
    for(i = 0; i < n; i++) {
        Up[i] = mesh->U[i] + K3[i] * dt;
    }
    diffusive_term(Up, n, h, d2U);
    convective_term(Up, n, h, dU);
    for(i = 0; i < n; i++) {
        K4[i] = DIFF_COEFF*d2U[i] - CONV_COEFF*dU[i];
    }

    // Compute Ui+1
    for(i = 0; i < n; i++) {
        mesh->U[i] += (K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) * dt/6.0;
    }

    // Free the memory
    free(dU);
    free(d2U);
    free(Up);
    free(K1);
    free(K2);
    free(K3);
    free(K4);
}

void compute_global_diagnostics(Mesh *mesh, double t, double *u, double *Q, double *E, double *R) {
    *Q = 0; *E = 0; *R = 0;
    for (int i = 0; i < mesh->n; i++) {
        *Q += mesh->U[i];
        *E += mesh->U[i] * mesh->U[i] / 2.0;
        *R += pow(mesh->U[i] - u[i], 2);
    }
    *Q *= mesh->h;
    *E *= mesh->h;
    *R = sqrt(mesh->h * (*R));
}

void compute_solution(double t, int n, int n_step, double *u, double *e) {
    *e = 1.0 / sqrt(8*M_PI * (pow(SIGMA_0, 2) + 4*DIFF_COEFF*t));
    for (int i = 0; i < n; i++) {
        u[(i+n_step) % n] = exp(-pow((i+n_step-n/2-1)*SPACE_DISCR - CONV_COEFF*t, 2)/(pow(SIGMA_0, 2) + 4*DIFF_COEFF*t)) / sqrt(M_PI * (pow(SIGMA_0, 2) + 4*DIFF_COEFF*t));
    }
}

Mesh *init_mesh() {
    Mesh *mesh = (Mesh*) malloc(sizeof(Mesh));
    mesh->n = N_DISCR;
    mesh->h = SPACE_DISCR;
    mesh->U = (double*) calloc(mesh->n, sizeof(double));

    // Initial condition u(x,0)
    for (int i = 0; i < mesh->n; i++) {
        mesh->U[i] = exp(-pow((i-mesh->n/2)*SPACE_DISCR, 2)/pow(SIGMA_0, 2)) / sqrt(pow(SIGMA_0,2)*M_PI);
    }   // NB: I chose Q = 1, arbitrarly.
    return mesh;
}

void free_mesh(Mesh *mesh) {
    free(mesh->U);
    free(mesh);
}
