#ifndef _DEBUG_H_
#define _DEBUG_H_

#include "mesh.h"
#include <math.h>

typedef struct {
    int n1, n2, n;
    double *val1, *val2;
    double dt;
} Tracker;

// Function declarations
void solenoidal_speed(double r, double theta, double *u, double *v);
Tracker *init_Tracker(double n1, double n2, double *val1, double *val2, double dt);
void free_Tracker(Tracker *t);
Tracker *track_mesh(Mesh *mesh, double dt);
void save_header(Tracker *t, double *x, double *y, FILE *fp);
void save_tracking_state(Tracker *t, double state, FILE *fp);

void set_speed_pressure(MACMesh *mesh, void (*f)(double, double, double*, double*));

#endif