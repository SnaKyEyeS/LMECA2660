#include "mesh.h"

MACMesh *init_mac_mesh(MappingType type) {
    MACMesh *mesh = (MACMesh *) malloc(sizeof(MACMesh));

    switch (type) {
        case AIRFOIL: {
            AirfoilMapping *mapping = init_airfoil_mapping();
            mesh->type = type;
            mesh->n1 = mapping->n_xi1;
            mesh->n2 = mapping->n_xi2 - 1;
            mesh->n_cell = mesh->n1 * mesh->n2;
            mesh->d1 = mapping->dxi1;
            mesh->d2 = mapping->dxi2;

            // Initialize the different meshes
            int ind;
            double xi1, xi2;
            mesh->w = init_mesh(mesh->n1+1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->w->n1; i++) {
                for (int j = 0; j < mesh->w->n2; j++) {
                    ind = i*mesh->w->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2;
                    airfoil_metrics(mapping, xi1, xi2, &mesh->w->x[ind], &mesh->w->y[ind], &mesh->w->h1[ind], &mesh->w->h2[ind],
                                    &mesh->w->dh1_d1[ind], &mesh->w->dh1_d2[ind], &mesh->w->dh2_d1[ind], &mesh->w->dh2_d2[ind],
                                    &mesh->w->ddh1_d1d2[ind], &mesh->w->ddh2_d1d2[ind], &mesh->w->theta[ind]);
                    mesh->w->val1[ind] = 0;
                    mesh->w->val2[ind] = 0;
                }
            }

            mesh->u = init_mesh(mesh->n1+1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->u->n1; i++) {
                for (int j = 0; j < mesh->u->n2; j++) {
                    ind = i*mesh->u->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2 + mapping->dxi2/2.0;
                    airfoil_metrics(mapping, xi1, xi2, &mesh->u->x[ind], &mesh->u->y[ind], &mesh->u->h1[ind], &mesh->u->h2[ind],
                                    &mesh->u->dh1_d1[ind], &mesh->u->dh1_d2[ind], &mesh->u->dh2_d1[ind], &mesh->u->dh2_d2[ind],
                                    &mesh->u->ddh1_d1d2[ind], &mesh->u->ddh2_d1d2[ind], &mesh->u->theta[ind]);
                    airfoil_init_velocity(mapping, U_INF, xi1, xi2, &mesh->u->val1[ind], NULL);
                    mesh->u->val2[ind] = 0;
                }
            }

            mesh->v = init_mesh(mesh->n1+1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->v->n1; i++) {
                for (int j = 0; j < mesh->v->n2; j++) {
                    ind = i*mesh->v->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1 + mapping->dxi1/2.0;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2;
                    airfoil_metrics(mapping, xi1, xi2, &mesh->v->x[ind], &mesh->v->y[ind], &mesh->v->h1[ind], &mesh->v->h2[ind],
                                    &mesh->v->dh1_d1[ind], &mesh->v->dh1_d2[ind], &mesh->v->dh2_d1[ind], &mesh->v->dh2_d2[ind],
                                    &mesh->v->ddh1_d1d2[ind], &mesh->v->ddh2_d1d2[ind], &mesh->v->theta[ind]);
                    airfoil_init_velocity(mapping, U_INF, xi1, xi2, NULL, &mesh->v->val1[ind]);
                    mesh->v->val2[ind] = 0;
                }
            }
            // The last column of v values is only to store a "ghost point.
            // Thus, we resize the mesh to its real value ! (mesh->v was initialized with n1 that was one bigger than required)
            mesh->v->n  -= mesh->v->n2;
            mesh->v->n1 -= 1;

            mesh->p = init_mesh(mesh->n1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->p->n1; i++) {
                for (int j = 0; j < mesh->p->n2; j++) {
                    ind = i*mesh->p->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1 + mapping->dxi1/2.0;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2 + mapping->dxi2/2.0;
                    airfoil_metrics(mapping, xi1, xi2, &mesh->p->x[ind], &mesh->p->y[ind], &mesh->p->h1[ind], &mesh->p->h2[ind],
                                    &mesh->p->dh1_d1[ind], &mesh->p->dh1_d2[ind], &mesh->p->dh2_d1[ind], &mesh->p->dh2_d2[ind],
                                    &mesh->p->ddh1_d1d2[ind], &mesh->p->ddh2_d1d2[ind], &mesh->p->theta[ind]);
                    airfoil_init_pressure(mapping, U_INF, xi1, xi2, &mesh->p->val1[ind]);
                    mesh->p->val2[ind] = 0;
                }
            }

            double dt_min_fourier = FOURIER_MAX * mapping->h_wall_normal *mapping->h_wall_normal / NU;
            double dt_min_CFL     = CFL_MAX * mapping->h_wall_normal / (2*U_INF);
            printf("dt fourier = %.10f\n", dt_min_fourier);
            printf("dt CFL     = %.10f\n", dt_min_CFL);
            mesh->dt = fmin(dt_min_fourier, dt_min_CFL);

            free(mapping);
            break;
        }

        case CYLINDER: {
            CylinderMapping *mapping = init_cylinder_mapping();
            mesh->type = type;
            //mesh->n1 = mapping->n_xi1 - 1;
            //mesh->n2 = mapping->n_xi2 - 1;
            mesh->n1 = mapping->n_xi1;
            mesh->n2 = mapping->n_xi2 - 1; // We don't go until 2PI
            mesh->n_cell = mesh->n1 * mesh->n2;
            mesh->d1 = mapping->dxi1;
            mesh->d2 = mapping->dxi2;

            // Initialize the different meshes
            int ind;
            double xi1, xi2;
            mesh->w = init_mesh(mesh->n1+1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->w->n1; i++) {
                for (int j = 0; j < mesh->w->n2; j++) {
                    ind = i*mesh->w->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2;
                    cylinder_metrics(mapping, xi1, xi2, &mesh->w->x[ind], &mesh->w->y[ind], &mesh->w->h1[ind], &mesh->w->h2[ind],
                                    &mesh->w->dh1_d1[ind], &mesh->w->dh1_d2[ind], &mesh->w->dh2_d1[ind], &mesh->w->dh2_d2[ind],
                                    &mesh->w->ddh1_d1d2[ind], &mesh->w->ddh2_d1d2[ind], &mesh->w->theta[ind]);
                    mesh->w->val1[ind] = 0;
                    mesh->w->val2[ind] = 0;
                }
            }

            mesh->u = init_mesh(mesh->n1+1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->u->n1; i++) {
                for (int j = 0; j < mesh->u->n2; j++) {
                    ind = i*mesh->u->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2 + mapping->dxi2/2.0;
                    cylinder_metrics(mapping, xi1, xi2, &mesh->u->x[ind], &mesh->u->y[ind], &mesh->u->h1[ind], &mesh->u->h2[ind],
                                    &mesh->u->dh1_d1[ind], &mesh->u->dh1_d2[ind], &mesh->u->dh2_d1[ind], &mesh->u->dh2_d2[ind],
                                    &mesh->u->ddh1_d1d2[ind], &mesh->u->ddh2_d1d2[ind], &mesh->u->theta[ind]);
                    cylinder_init_velocity(mapping, U_INF, xi1, xi2, &mesh->u->val1[ind], NULL);
                    mesh->u->val2[ind] = 0;
                }
            }

            mesh->v = init_mesh(mesh->n1+1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->v->n1; i++) {
                for (int j = 0; j < mesh->v->n2; j++) {
                    ind = i*mesh->v->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1 + mapping->dxi1/2.0;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2;
                    cylinder_metrics(mapping, xi1, xi2, &mesh->v->x[ind], &mesh->v->y[ind], &mesh->v->h1[ind], &mesh->v->h2[ind],
                                    &mesh->v->dh1_d1[ind], &mesh->v->dh1_d2[ind], &mesh->v->dh2_d1[ind], &mesh->v->dh2_d2[ind],
                                    &mesh->v->ddh1_d1d2[ind], &mesh->v->ddh2_d1d2[ind], &mesh->v->theta[ind]);
                    cylinder_init_velocity(mapping, U_INF, xi1, xi2, NULL, &mesh->v->val1[ind]);
                    mesh->v->val2[ind] = 0;

                }
            }
            // The last column of v values is only to store a "ghost point.
            // Thus, we resize the mesh to its real value ! (mesh->v was initialized with n1 that was one bigger than required)
            mesh->v->n  -= mesh->v->n2;
            mesh->v->n1 -= 1;

            mesh->p = init_mesh(mesh->n1, mesh->n2+1, mesh->d1, mesh->d2);
            for (int i = 0; i < mesh->p->n1; i++) {
                for (int j = 0; j < mesh->p->n2; j++) {
                    ind = i*mesh->p->n2 + j;
                    xi1 = mapping->xi1_lim[0] + i*mapping->dxi1 + mapping->dxi1/2.0;
                    xi2 = mapping->xi2_lim[0] + j*mapping->dxi2 + mapping->dxi2/2.0;
                    cylinder_metrics(mapping, xi1, xi2, &mesh->p->x[ind], &mesh->p->y[ind], &mesh->p->h1[ind], &mesh->p->h2[ind],
                                    &mesh->p->dh1_d1[ind], &mesh->p->dh1_d2[ind], &mesh->p->dh2_d1[ind], &mesh->p->dh2_d2[ind],
                                    &mesh->p->ddh1_d1d2[ind], &mesh->p->ddh2_d1d2[ind], &mesh->p->theta[ind]);
                    cylinder_init_pressure(mapping, U_INF, xi1, xi2, &mesh->p->val1[ind]);
                    mesh->p->val2[ind] = 0;
                }
            }

            double dt_min_fourier = FOURIER_MAX * mapping->h_wall_normal * mapping->h_wall_normal / NU;
            double dt_min_CFL     = CFL_MAX * mapping->h_wall_normal / (2*U_INF);
            printf("dt fourier = %.10f\n", dt_min_fourier);
            printf("dt CFL     = %.10f\n", dt_min_CFL);
            mesh->dt = fmin(dt_min_fourier, dt_min_CFL);

            free(mapping);
            break;
        }

        default:
            printf("Mesh type not recognized !\n");
            exit(0);
    }

    return mesh;
}


Mesh *init_mesh(int n1, int n2, double d1, double d2) {
    // Initialize the Mesh
    Mesh *mesh = (Mesh *) malloc(sizeof(Mesh));
    mesh->n = n1 * n2;
    mesh->n1 = n1;
    mesh->n2 = n2;
    mesh->d1 = d1;
    mesh->d2 = d2;

    // Allocate mesh data
    mesh->x = (double *) malloc(mesh->n*sizeof(double));
    mesh->y = (double *) malloc(mesh->n*sizeof(double));
    mesh->h1 = (double *) malloc(mesh->n*sizeof(double));
    mesh->h2 = (double *) malloc(mesh->n*sizeof(double));
    mesh->dh1_d1 = (double *) malloc(mesh->n*sizeof(double));
    mesh->dh1_d2 = (double *) malloc(mesh->n*sizeof(double));
    mesh->dh2_d1 = (double *) malloc(mesh->n*sizeof(double));
    mesh->dh2_d2 = (double *) malloc(mesh->n*sizeof(double));
    mesh->ddh1_d1d2 = (double *) malloc(mesh->n*sizeof(double));
    mesh->ddh2_d1d2 = (double *) malloc(mesh->n*sizeof(double));
    mesh->theta = (double *) malloc(mesh->n*sizeof(double));
    mesh->val1 = (double *) calloc(mesh->n, sizeof(double));
    mesh->val2 = (double *) calloc(mesh->n, sizeof(double));
}

void free_mac_mesh(MACMesh *mesh) {
    free_mesh(mesh->u);
    free_mesh(mesh->v);
    free_mesh(mesh->w);
    free_mesh(mesh->p);
}

void free_mesh(Mesh *mesh) {
    free(mesh->x);
    free(mesh->y);
    free(mesh->h1);
    free(mesh->h2);
    free(mesh->dh1_d1);
    free(mesh->dh1_d2);
    free(mesh->dh2_d1);
    free(mesh->dh2_d2);
    free(mesh->ddh1_d1d2);
    free(mesh->ddh2_d1d2);
    free(mesh->theta);
    free(mesh->val1);
    free(mesh->val2);
    free(mesh);
}

void save_array(double *arr, int size, FILE *fp) {
    for (int i = 0; i < size - 1; i++) {
        fprintf(fp, "%.10f, ", arr[i]);
    }
    fprintf(fp, "%.10f\n", arr[size-1]);
}

void save_mesh_header(Mesh *mesh, FILE *fp) {
    fprintf(fp, "{'n': %d, 'n1': %d, 'n2': %d, 'd1': %.10f, 'd2': %.10f}\n", mesh->n, mesh->n1, mesh->n2, mesh->d1, mesh->d2);
    save_array(mesh->x, mesh->n, fp);
    save_array(mesh->y, mesh->n, fp);
}

void save_mesh_state(Mesh *mesh, double state, FILE *fp) {
    fprintf(fp, "%.10f, ", state);
    save_array(mesh->val1, mesh->n, fp);
    fprintf(fp, "%.10f, ", state);
    save_array(mesh->val2, mesh->n, fp);
}

/*
 *  Saves the mesh to a file
 *  Column format:
 *  x, y, val1, val2
 */
void save_mesh(Mesh *mesh) {
    FILE *fp = fopen("../data/mesh.txt", "w+");
    for (int i = 0; i < mesh->n; i++) {
        fprintf(fp, "%.10f, %.10f, %.10f, %.10f\n", mesh->x[i], mesh->y[i], mesh->val1[i], mesh->val2[i]);
    }
    fclose(fp);
}
