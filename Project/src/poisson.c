#include "poisson.h"


/*Called by poisson_solver at each time step
  More than probably, you should need to add arguments to the prototype ...
  Modification to do :
      -Fill vector rhs*/
void computeRHS(Mesh *mesh, double *rhs, PetscInt rowStart, PetscInt rowEnd) {
    double r, theta, x, y;
    double r0 = 1;
    double r1 = 20+1;
    double R1R0 = 1/(r1-r0);
    double n=2;

    for(int i = rowStart; i < rowEnd ; i++) {
        x = mesh->x[i];
        y = mesh->y[i];
        r = sqrt(x*x + y*y);
        theta = atan2(y, x);

        rhs[i] = (1/r)*(M_PI*R1R0*sin(2*M_PI*(r-r0)*R1R0)*cos(n*theta))
						+ (2*pow(M_PI*R1R0,2)*cos(2*M_PI*(r-r0)*R1R0)*cos(n*theta))
						+ (1/(r*r))*(-pow(n,2)*pow(sin((r-r0)/(r1-r0)*M_PI),2)*cos(n*theta));
    }
}

/*To call at each time step after computation of U_star. This function solves the poisson equation
  and coM_PIes the solution of the equation into your vector Phi
  More than probably, you should need to add arguments to the prototype ...
  Modification to do :
      - Change the call to computeRHS as you have to modify its prototype too
      - Copy solution of the equation into your vector PHI*/
void poisson_solver(PoissonData *data, Mesh *mesh)
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
        mesh->val2[i] = sol[i];
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
void computeLaplacianMatrix(Mesh *mesh, Mat A, int rowStart, int rowEnd) {
    double d1 = mesh->d1;
    double d2 = mesh->d2;
    double h1, h2, dh1_d1, dh1_d2, dh2_d1, dh2_d2, a, b;
    double phi, phi_i_plus, phi_i_minus, phi_j_plus, phi_j_minus;
    int i, j;

    // for (int i = rowStart; i < rowEnd; i++) {
    //     MatSetValue(A, i, i , 1.0, INSERT_VALUES);  // Main diagonal
    // }

    for (int r = rowStart; r < rowEnd; r++){
        i = r / mesh->n2;
        j = r % mesh->n2;

        h1 = mesh->h1[r];
        h2 = mesh->h2[r];
        dh1_d1 = mesh->dh1_d1[r];
        dh1_d2 = mesh->dh1_d2[r];
        dh2_d1 = mesh->dh2_d1[r];
        dh2_d2 = mesh->dh2_d2[r];
        a = (dh2_d1*h1 - h2*dh1_d1)/(h1*h1*h1*h2);
        b = (dh1_d2*h2 - h1*dh2_d2)/(h2*h2*h2*h1);

        phi_i_plus = 1/(d1*d1*h1*h1) + a/(2*d1);
        phi_i_minus = 1/(d1*d1*h1*h1) - a/(2*d1);
        phi_j_plus = 1/(d2*d2*h2*h2) + b/(2*d2);
        phi_j_minus = 1/(d2*d2*h2*h2) - b/(2*d2);

        if (i == 0) {
            phi = -1/(d1*d1*h1*h1) -2/(d2*d2*h2*h2) + a/(2*d1);
        } else if (i == mesh->n1-1) {
            phi = -1/(d1*d1*h1*h1) -2/(d2*d2*h2*h2) - a/(2*d1);
        } else {
            phi = -2/(d1*d1*h1*h1) -2/(d2*d2*h2*h2);
        }

        // Main diagonal
        MatSetValue(A, r, r, phi, INSERT_VALUES);

        // Right & left borders
        if (i > 0) {
            MatSetValue(A, r, r-mesh->n2, phi_i_minus, INSERT_VALUES);
        }
        if (i < mesh->n1-1) {
            MatSetValue(A, r, r+mesh->n2, phi_i_plus, INSERT_VALUES);
        }

        // Upper & lower borders (those are periodic)
        if (j > 0) {
            MatSetValue(A, r, r-1, phi_j_minus, INSERT_VALUES);
        } else {
            MatSetValue(A, r, r-1+mesh->n2, phi_j_minus, INSERT_VALUES);
        }
        if (j < mesh->n2-1) {
            MatSetValue(A, r, r+1, phi_j_plus, INSERT_VALUES);
        } else {
            MatSetValue(A, r, r+1-mesh->n2, phi_j_plus, INSERT_VALUES);
        }
    }
}

/*To call during the initialization of your solver, before the begin of the time loop
  Maybe you should need to add an argument to specify the number of unknows
  Modification to do in this function :
      -Specify the number of unknows
      -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(PoissonData* data, Mesh *mesh) {
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
    if (0) {
        PetscViewer viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, "../data/test.m", &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView(data->A,viewer);
        PetscViewerPopFormat(viewer);
        PetscViewerDestroy(&viewer);
    }


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
