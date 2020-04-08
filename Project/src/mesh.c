#include "mesh.h"


OrthogonalMesh *init_mesh(MappingType type) {
    // Initialize the Mesh
    OrthogonalMesh *mesh = (OrthogonalMesh *) malloc(sizeof(OrthogonalMesh));

    // Initialize the mesh
    switch (type) {
        case AIRFOIL: {     // Airfoil - Joukowsky mapping
            AirfoilMapping *mapping = init_airfoil_mapping();
            mesh->type = type;
            mesh->n_xi1 = mapping->n_xi1;
            mesh->n_xi2 = mapping->n_xi2;
            mesh->n = mesh->n_xi1 * mesh->n_xi2;
            mesh->dxi1 = mapping->dxi1;
            mesh->dxi2 = mapping->dxi2;

            // Allocate mesh data
            mesh->data = (double *) malloc(14*mesh->n*sizeof(double));
            mesh->x = &mesh->data[0*mesh->n];
            mesh->y = &mesh->data[1*mesh->n];
            mesh->h1 = &mesh->data[2*mesh->n];
            mesh->h2 = &mesh->data[3*mesh->n];
            mesh->dh1_d1 = &mesh->data[4*mesh->n];
            mesh->dh1_d2 = &mesh->data[5*mesh->n];
            mesh->dh2_d1 = &mesh->data[6*mesh->n];
            mesh->dh2_d2 = &mesh->data[7*mesh->n];
            mesh->ddh1_d1d2 = &mesh->data[8*mesh->n];
            mesh->ddh2_d1d2 = &mesh->data[9*mesh->n];
            mesh->theta = &mesh->data[10*mesh->n];
            mesh->u_n = &mesh->data[11*mesh->n];
            mesh->u_theta = &mesh->data[12*mesh->n];
            mesh->p = &mesh->data[13*mesh->n];

            // Initialize mesh data
            int ind;
            double xi1, xi2;
            for (int i = 0; i < mesh->n_xi1; i++) {
                for (int j = 0; j < mesh->n_xi2; j++) {
                    ind = i*mesh->n_xi2 + j;
                    xi1 = mapping->xi1_lim[0] + (mapping->xi1_lim[1]-mapping->xi1_lim[0])*i*mapping->dxi1;
                    xi2 = mapping->xi2_lim[0] + (mapping->xi1_lim[1]-mapping->xi1_lim[0])*j*mapping->dxi2;
                    airfoil_metrics(mapping, xi1, xi2, &mesh->x[ind], &mesh->y[ind], &mesh->h1[ind], &mesh->h2[ind],
                                    &mesh->dh1_d1[ind], &mesh->dh1_d2[ind], &mesh->dh2_d1[ind], &mesh->dh2_d2[ind],
                                    &mesh->ddh1_d1d2[ind], &mesh->ddh2_d1d2[ind], &mesh->theta[ind]);
                    airfoil_init_velocity(mapping, U_INF, xi1, xi2, &mesh->u_n[ind], &mesh->u_theta[ind]);
                    airfoil_init_pressure(mapping, U_INF, xi1, xi2, &mesh->p[ind]);
                }
            }

            free(mapping);
            break; }

        case CYLINDER: {    // Cylinder
            CylinderMapping *mapping = init_cylinder_mapping();
            mesh->type = type;
            mesh->n_xi1 = mapping->n_xi1;
            mesh->n_xi2 = mapping->n_xi2;
            mesh->n = mesh->n_xi1 * mesh->n_xi2;
            mesh->dxi1 = mapping->dxi1;
            mesh->dxi2 = mapping->dxi2;

            // Allocate mesh data
            mesh->data = (double *) malloc(14*mesh->n*sizeof(double));
            mesh->x = &mesh->data[0*mesh->n];
            mesh->y = &mesh->data[1*mesh->n];
            mesh->h1 = &mesh->data[2*mesh->n];
            mesh->h2 = &mesh->data[3*mesh->n];
            mesh->dh1_d1 = &mesh->data[4*mesh->n];
            mesh->dh1_d2 = &mesh->data[5*mesh->n];
            mesh->dh2_d1 = &mesh->data[6*mesh->n];
            mesh->dh2_d2 = &mesh->data[7*mesh->n];
            mesh->ddh1_d1d2 = &mesh->data[8*mesh->n];
            mesh->ddh2_d1d2 = &mesh->data[9*mesh->n];
            mesh->theta = &mesh->data[10*mesh->n];
            mesh->u_n = &mesh->data[11*mesh->n];
            mesh->u_theta = &mesh->data[12*mesh->n];
            mesh->p = &mesh->data[13*mesh->n];

            // Initialize mesh data
            int ind;
            double xi1, xi2;
            for (int i = 0; i < mesh->n_xi1; i++) {
                for (int j = 0; j < mesh->n_xi2; j++) {
                    ind = i*mesh->n_xi2 + j;
                    xi1 = mapping->xi1_lim[0] + (mapping->xi1_lim[1]-mapping->xi1_lim[0])*i*mapping->dxi1;
                    xi2 = mapping->xi2_lim[0] + (mapping->xi1_lim[1]-mapping->xi1_lim[0])*j*mapping->dxi2;
                    cylinder_metrics(mapping, xi1, xi2, &mesh->x[ind], &mesh->y[ind], &mesh->h1[ind], &mesh->h2[ind],
                                    &mesh->dh1_d1[ind], &mesh->dh1_d2[ind], &mesh->dh2_d1[ind], &mesh->dh2_d2[ind],
                                    &mesh->ddh1_d1d2[ind], &mesh->ddh2_d1d2[ind], &mesh->theta[ind]);
                    cylinder_init_velocity(mapping, U_INF, xi1, xi2, &mesh->u_n[ind], &mesh->u_theta[ind]);
                    cylinder_init_pressure(mapping, U_INF, xi1, xi2, &mesh->p[ind]);
                }
            }

            free(mapping);
            break; }

        default:
            printf("Mesh type not recognized !\n");
            exit(0);
    }

    return mesh;
}

void free_mesh(OrthogonalMesh *mesh) {
    free(mesh->data);
    free(mesh);
}

/*
 *  Saves the mesh to a file
 *  Column format:
 *  x, y, u_n, u_theta, p
 */
void save_mesh(OrthogonalMesh *mesh) {
    FILE *fp = fopen("../data/mesh.txt", "w+");
    for (int i = 0; i < mesh->n; i++) {
        fprintf(fp, "%f, %f, %f, %f, %f\n", mesh->x[i], mesh->y[i], mesh->u_n[i], mesh->u_theta[i], mesh->p[i]);
    }
    fclose(fp);
}
