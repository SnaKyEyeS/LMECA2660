#include "poisson.h"


/*Called by poisson_solver at each time step
  More than probably, you should need to add arguments to the prototype ...
  Modification to do :
      -Fill vector rhs*/
void computeRHS(OrthogonalMesh *mesh, double *rhs, PetscInt rowStart, PetscInt rowEnd) {
    double r, theta, x, y;
    for(int i = rowStart; i < rowEnd ; i++) {
        x = mesh->x[i];
        y = mesh->y[i];
        r = sqrt(x*x + y*y);
        if (r < 0.1) {
            printf("caca boudin prout pipi caca\n");
            r = 0.1;
        }
        theta = atan2(y, x);

        double s = sin(M_PI*(0.1-r)/10);
        double c = cos(2*theta);
        rhs[i] = s*c*M_PI*M_PI / (10*10) + cos(M_PI*(0.1-r)/10)*c*M_PI / (10*r) + 4*s*c / (r*r);
    }
}

/*To call at each time step after computation of U_star. This function solves the poisson equation
  and coM_PIes the solution of the equation into your vector Phi
  More than probably, you should need to add arguments to the prototype ...
  Modification to do :
      - Change the call to computeRHS as you have to modify its prototype too
      - Copy solution of the equation into your vector PHI*/
void poisson_solver(PoissonData *data, OrthogonalMesh *mesh)
{

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = data->sles;
    Vec b = data->b;
    Vec x = data->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(mesh, rhs, rowStart, rowEnd); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    for(int i = rowStart; i < rowEnd; i++){
        mesh->p[i] = sol[i];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.
  In its current state, it inserts unity on the main diagonal.
  More than probably, you should need to add arguments to the prototype ... .
  Modification to do in this function :
      -Insert the correct factor in matrix A*/
void computeLaplacianMatrix(OrthogonalMesh *mesh, Mat A, int rowStart, int rowEnd) {
    double dx_sq = mesh->dxi1*mesh->dxi1;
    double dy_sq = mesh->dxi2*mesh->dxi2;

    for (int i = rowStart; i < rowEnd; i++) {
        MatSetValue(A, i, i , 1.0, INSERT_VALUES);  // Main diagonal
    }
    // for (int i = rowStart; i < rowEnd; i++){
    //     if (i < mesh->n_xi1)            { continue; }   // Lower border
    //     if (i > mesh->n - mesh->n_xi1)  { continue; }   // Upper border
    //     // if ((i % mesh->n_xi1) == 0)     {
    //     //
    //     // }   // Left  border
    //     // if (((i+1) % mesh->n_xi1) == 0) { continue; }   // Right border
    //
    //     // Main diagonal
    //     MatSetValue(A, i, i , -2*(1/dx_sq + 1/dy_sq), INSERT_VALUES);
    //
    //     // Left & Right elements
    //     if (i < mesh->n_xi1) {                          // Left border
    //         // MatSetValue(A, i, i+1, 1/dx_sq, INSERT_VALUES);
    //         // MatSetValue(A, i, i+mesh->n_xi1-1, 1/dx_sq, INSERT_VALUES);
    //
    //     } else if (i > mesh->n - mesh->n_xi1) {         // Right border
    //         // MatSetValue(A, i, i-mesh->n_xi1+1, 1/dx_sq, INSERT_VALUES);
    //         // MatSetValue(A, i, i-1, 1/dx_sq, INSERT_VALUES);
    //
    //     } else {
    //         MatSetValue(A, i, i+1, 1/dx_sq, INSERT_VALUES);
    //         MatSetValue(A, i, i-1, 1/dx_sq, INSERT_VALUES);
    //     }
    //
    //     // Up & Down elements
    //     MatSetValue(A, i, i+mesh->n_xi2, 1/dy_sq, INSERT_VALUES);
    //     MatSetValue(A, i, i-mesh->n_xi2, 1/dy_sq, INSERT_VALUES);
    // }
}

/*To call during the initialization of your solver, before the begin of the time loop
  Maybe you should need to add an argument to specify the number of unknows
  Modification to do in this function :
      -Specify the number of unknows
      -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(PoissonData* data, OrthogonalMesh *mesh) {
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	int nphi = mesh->n; /*WRITE HERE THE NUMBER OF UNKNOWS*/

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A, 5, NULL); // /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(mesh, data->A, rowStart, rowEnd);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Singular matrix */
    MatNullSpace nullspace;
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
    MatSetNullSpace(data->A, nullspace);
    MatSetTransposeNullSpace(data->A, nullspace);
    MatNullSpaceDestroy(&nullspace);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    PC prec;
    KSPGetPC(data->sles, &prec);
    KSPSetFromOptions(data->sles); // specifier le solveur exact a partir de la ligne commande mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-10, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,5);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Matrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation
  Modification to do : nothing */
void free_poisson_solver(PoissonData* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
