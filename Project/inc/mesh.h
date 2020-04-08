#ifndef _MESH_H_
#define _MESH_H_

#include <stdlib.h>
#include <stdio.h>
#include "airfoil.h"
#include "cylinder.h"

#define U_INF 1.0

typedef enum {AIRFOIL, CYLINDER} MappingType;
typedef struct {
    // Mapping
    MappingType type;
    void *mapping;

    // Number of unknowns
    int n, n_xi1, n_xi2;

    // Mesh data
    double *data;
    double *x, *y;
    double *h1, *h2;
    double *dh1_d1, *dh1_d2, *dh2_d1, *dh2_d2;
    double *ddh1_d1d2, *ddh2_d1d2;
    double *theta;
    double *u_n, *u_theta;
    double *p;

} OrthogonalMesh;


OrthogonalMesh *init_mesh(MappingType type);
void            free_mesh(OrthogonalMesh *mesh);
void            save_mesh(OrthogonalMesh *mesh);


#endif
