#ifndef _MESH_H_
#define _MESH_H_

#include <stdlib.h>
#include <stdio.h>
#include "airfoil.h"
#include "cylinder.h"

#define N_MESH 4

#define REYNOLDS    550
#define FOURIER     0.1
#define CFL         0.7

#define NU          (1e-5)      // UNUSED VARIABLE -> defined for the save_header
                                // in debug.c --> to change eventually !

typedef enum {AIRFOIL, CYLINDER} MappingType;



typedef struct {
    // Number of unknowns & discretization
    int n;          // Number of nodes
    int n1, n2;     // Number of nodes in each direction
    double d1, d2;  // Step between each node

    // Mesh data
    double *x, *y;
    double *h1, *h2;
    double *dh1_d1, *dh1_d2, *dh2_d1, *dh2_d2;
    double *ddh1_d1d2, *ddh2_d1d2;
    double *theta;
    double *val1, *val2;        // Can hold several values if needed.

} Mesh;

typedef struct {
    MappingType type;   // Mapping type
    int n_cell;         // Total number of cells
    int n1, n2;         // Number of cell in the (1) and (2) direction
    double d1, d2;      // Size of the cells

    double dt;          // Timestep for the numerical integration scheme
    double Lc;          // Mesh's caracteristic length
    double Uinf;
    double nu;

    Mesh *w;
    Mesh *u;
    Mesh *v;
    Mesh *p;
} MACMesh;

MACMesh    *init_mac_mesh(MappingType type);
void        free_mac_mesh(MACMesh *mesh);

Mesh   *init_mesh(int n1, int n2, double d1, double d2);
void    free_mesh(Mesh *mesh);
void    save_array(double *arr, int size, FILE *fp);
void    save_mesh_header(Mesh *mesh, FILE *fp);
void    save_mesh_state(Mesh *mesh, double state, FILE *fp);
void    save_mesh(Mesh *mesh);


#endif
