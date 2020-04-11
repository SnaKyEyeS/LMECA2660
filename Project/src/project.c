#include "poisson.h"
#include "airfoil.h"
#include "cylinder.h"
#include "mesh.h"
#include "stdio.h"

/*Example : discretization of \nabla \cdot \mathbf{u}. Do not compile this code, it will not; that's a example.
  In this example, I have :
     - Data : Structure defined by myself (you should do the same) that contains the important data for
              solver, such as 1) the number of point on the xi1 and xi2 axis for the P, Un and Utheta field,
              2) the resolution dxi1 and dxi2 and 3) tables refecencing form factors and their derivatives.
     - Un : a pointer referencing to a table containing the normal component of the velocity field.
     - Utheta : a pointer referencing to a table containing the tangential component of the velocity field.
     - div_u : a pointer referencing to a table that will be filled with the divergence of v = (Un, Utheta) */
// void divergence(Data *data, double ** Un, double ** Utheta, double ** div_u){

//     int nxi1_P = data->nxi1_P;
//     int nxi2_P = data->nxi2_P;

//     int nxi1_Ua = data->nxi1_Ua;
//     int nxi2_Ua = data->nxi2_Ua;

//     int nxi1_Ub = data->nxi1_Ub;
//     int nxi2_Ub = data->nxi2_Ub;

//     double ** h2_p = data->h2_p;
//     double ** h1_p = data->h1_p;
//     double ** h2_un = data->h2_un;
//     double ** h1_utheta = data->h1_utheta;

//     double dxi1 = data->dxi1;
//     double dxi2 = data->dxi2;

//     int i,j;
//     for(j=1;j<(nxi2_P-1);j++){
//         for(i=1;i<(nxi1_P-1);i++){
//             div_u[j][i] = (1.0/(h2_p[j][i]*h1_p[j][i]))*( (h2_un[j][i]*Un[j][i] - h2_un[j][i-1]*Un[j][i-1])/dxi1 \
//                           + (h1_utheta[j][i]*Utheta[j][i] - h1_utheta[j-1][i]*Utheta[j-1][i])/dxi2 );
//         }
//     }
// }


int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);

    // Initialize Mesh
    MACMesh *mesh = init_mac_mesh(CYLINDER);

    // Initialize Poisson solver
    PoissonData *poisson = (PoissonData *) malloc(sizeof(PoissonData));
    initialize_poisson_solver(poisson, mesh->p);


    // Compute Poisson Solution
    poisson_solver(poisson, mesh->p);
    save_mesh(mesh->p);

    // Free memory
    free_mac_mesh(mesh);
    free_poisson_solver(poisson);
    PetscFinalize();

}
