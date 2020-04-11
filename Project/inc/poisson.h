#ifndef _POISSON_H_
#define _POISSON_H_

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include "mesh.h"

//Structure storing petsc vectors
typedef struct {
	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} PoissonData;

PetscErrorCode initialize_poisson_solver(PoissonData* data, Mesh *mesh);
void poisson_solver(PoissonData *data, Mesh *mesh);
void free_poisson_solver(PoissonData* data);

#endif
