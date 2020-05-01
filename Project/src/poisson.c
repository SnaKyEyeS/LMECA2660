#include "poisson.h"
#include "solver.h"


/*Called by poisson_solver at each time step
  More than probably, you should need to add arguments to the prototype ...
  Modification to do :
      -Fill vector rhs*/
void computeRHS(MACMesh *mesh, double *rhs, PetscInt rowStart, PetscInt rowEnd) {
    double r, theta, x, y;
    double r0 = 1.0/100.0;
    double r1 = 1 + r0;
    double n=2;
    double fact;

    for(int ind = rowStart; ind < rowEnd ; ind++) {

        if (ind == 0) {
            rhs[ind] =  0.0;
            continue;
        }

        x = mesh->p->x[ind];
        y = mesh->p->y[ind];
        r = sqrt(x*x + y*y);
        theta = atan2(y, x);

        fact = M_PI * (r0 - r) / (r1 - r0);
        rhs[ind] = 2*cos(2*fact)*cos(2*theta)*M_PI*M_PI / pow(r1-r0, 2)
                    - 2*cos(fact)*sin(fact)*cos(2*theta)*M_PI / ((r1-r0)*r)
                    - 4*pow(sin(fact), 2)*cos(2*theta) / pow(r, 2);
    }
}

/*To call at each time step after computation of U_star. This function solves the poisson equation
  and coM_PIes the solution of the equation into your vector Phi
  More than probably, you should need to add arguments to the prototype ...
  Modification to do :
      - Change the call to computeRHS as you have to modify its prototype too
      - Copy solution of the equation into your vector PHI*/
void poisson_solver(PoissonData *data, MACMesh *mesh)
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
    // computeRHS(mesh, rhs, rowStart, rowEnd); // This is to test the poisson solver !
    compute_rhs(mesh, rhs, mesh->dt);
    VecRestoreArray(b, &rhs);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, " \tPoisson solved in %d iterations.\n", its);

    VecGetArray(x, &sol);

    for(int i = rowStart; i < rowEnd; i++){
        mesh->p->val2[i] = sol[i];
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
void computeLaplacianMatrix(MACMesh *mesh, Mat A, int rowStart, int rowEnd) {
    double d1 = mesh->p->d1;
    double d2 = mesh->p->d2;
    double h1, h2, dh1_d1, dh1_d2, dh2_d1, dh2_d2, a, b;
    double phi, phi_right, phi_left, phi_up, phi_bottom;
    double h1h2;
    double form_right, form_left, form_up, form_bottom;

    int ind;
    int ind_phi_left, ind_phi_right, ind_phi_bottom, ind_phi_up;
    int ind_u_left, ind_u_right;
    int ind_v_bottom, ind_v_up;

    double i_den, j_den;

    for (int i = 0; i < mesh->p->n1; i++) {
        for (int j = 0; j < mesh->p->n2; j++) {

            ind = i * mesh->p->n2 + j;
            if (ind == 0) {
                MatSetValue(A, ind, ind, 1.0, INSERT_VALUES);
                continue;
            }

            ind_phi_left    = index(i, j, mesh->p->n2, -1, 0);
            ind_phi_right   = index(i, j, mesh->p->n2, 1, 0);
            ind_phi_bottom  = index(i, j, mesh->p->n2, 0, -1);
            ind_phi_up      = index(i, j, mesh->p->n2, 0, 1);

            ind_u_left      = index(i, j, mesh->u->n2, 0, 0);
            ind_u_right     = index(i, j, mesh->u->n2, 1, 0);

            ind_v_bottom    = index(i, j, mesh->v->n2, 0, 0);
            ind_v_up        = index(i, j, mesh->v->n2, 0, 1);

            h1h2 = mesh->p->h1[ind] * mesh->p->h2[ind];
            i_den = 1 / (h1h2*d1*d1);
            j_den = 1 / (h1h2*d2*d2);

            form_right  = mesh->u->h2[ind_u_right]  / mesh->u->h1[ind_u_right];
            form_left   = mesh->u->h2[ind_u_left]   / mesh->u->h1[ind_u_left];
            form_bottom = mesh->v->h1[ind_v_bottom] / mesh->v->h2[ind_v_bottom];
            form_up     = mesh->v->h1[ind_v_up]     / mesh->v->h2[ind_v_up];

            phi_right   = form_right    * i_den;
            phi_left    = form_left     * i_den;
            phi_up      = form_up       * j_den;
            phi_bottom  = form_bottom   * j_den;


            if (i == 0) {
                phi =   -(form_right    + 0.0           ) * i_den
                        -(form_up       + form_bottom   ) * j_den;
            } else if (i == mesh->p->n1-1) {
                phi =   -(0.0           + form_left     ) * i_den
                        -(form_up       + form_bottom   ) * j_den;
            } else {
                phi =   -(form_right    + form_left     ) * i_den
                        -(form_up       + form_bottom   ) * j_den;
            }

            // Main diagonal
            MatSetValue(A, ind, ind, phi, INSERT_VALUES);

            // Right & left borders
            if (i > 0) {
                MatSetValue(A, ind, ind_phi_left, phi_left, INSERT_VALUES);
            }
            if (i < mesh->p->n1-1) {
                MatSetValue(A, ind, ind_phi_right, phi_right, INSERT_VALUES);
            }

            // Upper & lower borders (those are periodic)
            MatSetValue(A, ind, ind_phi_bottom, phi_bottom, INSERT_VALUES);
            MatSetValue(A, ind, ind_phi_up, phi_up, INSERT_VALUES);
        }
    }
}

/*To call during the initialization of your solver, before the begin of the time loop
  Maybe you should need to add an argument to specify the number of unknows
  Modification to do in this function :
      -Specify the number of unknows
      -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(PoissonData* data, MACMesh *mesh) {
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	int nphi = mesh->p->n; /*WRITE HERE THE NUMBER OF UNKNOWS*/

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
    MatSeqAIJSetPreallocation(data->A, 7, NULL); // /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(mesh, data->A, rowStart, rowEnd);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);


    if(0) {
        PetscViewer viewer;
    	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"../data/test.output",&viewer);
    	MatView(data->A,viewer);
    	PetscViewerDestroy(&viewer);
    }

    /* Singular matrix */
    // MatNullSpace nullspace;
    // MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
    // MatSetNullSpace(data->A, nullspace);
    // MatSetTransposeNullSpace(data->A, nullspace);
    // MatNullSpaceDestroy(&nullspace);

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
