#include <stdio.h>
#include <math.h>
#include <string.h>
#include "solver.h"


void save_data(char *filename, double *data, int i_end, int j_end) {
    char prefix[50] = "../matlab/data/";
    strcat(prefix, filename);
    FILE *fp = fopen(prefix, "w+");
    for (int i = 0; i < i_end; i++) {
        for (int j = 0; j < j_end; j++) {
            fprintf(fp, "%f, ", data[i*j_end+j]);
        }
        fprintf(fp, "\n");
    }
}

int main() {
    // Initialize stuff
    int n_step = T_END/TIME_DISCR;
    double t = 0;
    double dt = TIME_DISCR;
    Mesh *mesh = init_mesh();

    double *Q = (double*) calloc(n_step*3, sizeof(double));
    double *E = &Q[n_step];
    double *R = &Q[n_step*2];

    double *uh = (double*) calloc(mesh->n*n_step, sizeof(double));
    double *u = (double*) calloc(mesh->n*n_step, sizeof(double));
    double *e = (double*) calloc(n_step, sizeof(double));
    double e_0 = 0;
    compute_solution(t, mesh->n, 0, u, &e_0);

    for (int i = 0; i < n_step; i++) {
        // Do one RK4 iteration
        runge_kutta_4(mesh, dt);
        for (int j = 0; j < mesh->n; j++) {
            uh[i*mesh->n+j] = mesh->U[j];
        }

        // Compute exact solution
        compute_solution(t, mesh->n, i, &u[i*mesh->n], &e[i]);
        e[i] /= e_0;    // Normalize the enery

        // Compute global diagnostics
        compute_global_diagnostics(mesh, t, &u[i*mesh->n], &Q[i], &E[i], &R[i]);
        E[i] /= e_0;    // Normalize the enery (NB: Q is already normalized)
                        // since Q = 1 (arbitrarly fixed)
        R[i] /= sqrt(e_0);

        // Increment t
        t += dt;
    }

    // Save data
    save_data("uh.txt", uh, n_step, mesh->n);
    save_data("u.txt", u, n_step, mesh->n);
    save_data("global_diagnostics.txt", Q, 3, n_step);

    // Save data (run by run)
    // save_data("convection-diffusion/uh_125_i6.txt", uh, n_step, mesh->n);
    // save_data("convection-diffusion/u_125.txt", uh, n_step, mesh->n);
    // save_data("convection-diffusion/glob_125_i6.txt", Q, 3, n_step);

    // Save data for order of convergence
    // FILE *fp = fopen("../matlab/data/I6.txt", "a");
    // fprintf(fp, "%.32f, %.32f\n", mesh->h/SIGMA_0, R[n_step-1]);

    // Free the memory !
    free(Q);
    free(u);
    free_mesh(mesh);
    free(uh);
    free(e);

    return 0;
}
