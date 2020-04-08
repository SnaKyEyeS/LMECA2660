#ifndef _AIRFOIL_H_
#define _AIRFOIL_H_

#include <iostream>
#include <complex>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <limits>

typedef struct {

	//IMPORTANT FOR STUDENTS : number of cell without the ghosts in the normal and tangential direction.
	int n_xi1;
	int n_xi2;

	//IMPORTANT FOR STUDENTS : discretization in the computational domain.
	double dxi1;
	double dxi2;

	//IMPORTANT FOR STUDENTS : chord length, useful to make your fields dimensionless...
	double c; //Chord length.

	//NOT IMPORTANT FOR STUDENTS
	double fact1; //Factor alpha
	double fact2; //Factor beta
	double L_carac; //Diameter of the cylinder before the Joukowski tranformation.
	double Str_normal; //Facteur de stretching in the normal direction before the Joukowski tranformation.
	double h_wall_normal; //Fist size cell of the cylinder (before the Joukowski tranformation).
	double epsOa; //ratio between espsilon and the radius of the cylinder (before Joukowski tranformation)
	double reg; //Regularization factor.
	double bj; //b factor in the second order equation.
	double AoA; //AoA of the airfoil.
	double delta; //Camber parameter.
	double L; //Length of the domain in the xi_1 direction.
	double H; //Length of the domain in the xi_2 direction.
	double xi1_lim[2]; //Limit values for xi_1 in the computational domain.
	double xi2_lim[2]; //Limit values for x2_1 in the computational domain.


} AirfoilMapping;


// Function declarations
AirfoilMapping *init_airfoil_mapping();

/*PRE : - Mapping must be initialized with function airfoil_init_mapping before calling this function
		- xi1 and xi2 are mesh point coordinates the computational domain.
  POST : Based on the Joukowski transformation, this function computes :
  			- x, y : mesh point coordinates in the physical domain
			- h1, h2: form factors
			- dh1_dxi1, dh1_dxi2 : derivatives of the form factors h1 with respect to xi1 and xi2
			- d2h1_dxi1dxi2, d2h2_dxi1dxi2 : cross derivatives of h1 and h2
			- theta_mesh : local orientation of the mesh in the physical domain
		 If a information is not needed, replace its argument by NULL
		    Ex : theta_mesh is not needed => theta_mesh = NULL*/
void airfoil_metrics(AirfoilMapping* mapping, double xi1, double xi2, double * x, double * y, \
							  double * h1, double * h2, double * dh1_dxi1, double * dh1_dxi2,\
							  double * dh2_dxi1, double * dh2_dxi2, \
							  double * d2h1_dxi1dxi2, double * d2h2_dxi1dxi2, \
							  double * theta_mesh);

/*PRE : - Mapping must be initialized with function airfoil_init_mapping before calling this function
		- xi1 and xi2 are mesh point coordinates the computational domain.
		- U_inf is the upstream velocity
  POST : Based on the Joukowski transformation, this function computes :
  			- x, y : mesh point coordinates in the physical domain
			- u_n, u_theta : the normal and tangential component velocity of a potential flow.
		 If a information is not needed, replace its argument by NULL
		    Ex : u_n is not needed => u_n = NULL*/
void airfoil_init_velocity(AirfoilMapping* mapping, double U_inf, double xi1, double xi2, double * u_n, double *u_theta);

/*PRE : - Mapping must be initialized with function airfoil_init_mapping before calling this function
		- xi1 and xi2 are mesh point coordinates the computational domain.
		- U_inf is the upstream velocity
  POST : Based on the Joukowski transformation, this function computes :
  			- x, y : mesh point coordinates in the physical domain
			- u_n, u_theta : the normal and tangential component velocity of a potential flow*/
void airfoil_init_pressure(AirfoilMapping* mapping, double U_inf, double xi1, double xi2, double * P);

#endif
