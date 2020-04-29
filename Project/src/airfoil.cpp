#include "airfoil.h"

AirfoilMapping *init_airfoil_mapping() {
	AirfoilMapping *mapping = (AirfoilMapping*) malloc(sizeof(AirfoilMapping));

	mapping->Str_normal = 1.015;
	mapping->h_wall_normal = 2.0e-4;
	mapping->L_carac = 0.025;
	double n_D = 50.0;

	mapping->fact2 = mapping->h_wall_normal/(mapping->Str_normal-1);
	mapping->fact1 = log(((n_D*mapping->L_carac)/mapping->fact2) + 1);
	mapping->dxi1 = log(mapping->Str_normal)/mapping->fact1;
	mapping->n_xi1 = ceil(1.0/mapping->dxi1);

	mapping->epsOa = 0.15; //To set the width
	mapping->delta = M_PI/4.0; //To set the camber
	mapping->reg = 1.05; //To make the airfoil regularized at the leading edge.
	mapping->AoA = 20;
	double a = (mapping->L_carac/2.0);
	double eps = (mapping->epsOa)*a;
	double Delta = 4.0*eps*eps*cos(mapping->delta)*cos(mapping->delta) - 4.0*(eps*eps - a*a);
	mapping->bj = (-2.0*eps*cos(mapping->delta) + sqrt(Delta) )/2.0;

	mapping->n_xi2 = 360;

	mapping->L = mapping->n_xi1*mapping->dxi1;
	mapping->H = 2*M_PI;
	mapping->xi1_lim[0] = 0.0;
	mapping->xi2_lim[0] = 0.0;

	mapping->xi1_lim[1] = mapping->xi1_lim[0] + mapping->L;
	mapping->xi2_lim[1] = mapping->xi2_lim[0] + mapping->H;
	mapping->dxi2 = mapping->H/mapping->n_xi2;

	mapping->Lc = mapping->reg*4.0*mapping->bj;
	return mapping;
}




void airfoil_metrics(AirfoilMapping* mapping, double xi1, double xi2, double * x, double * y, \
							  double * h1, double * h2, double * dh1_dxi1, double * dh1_dxi2,\
							  double * dh2_dxi1, double * dh2_dxi2, \
							  double * d2h1_dxi1dxi2, double * d2h2_dxi1dxi2, \
							  double * theta_mesh){

	double fact1 = mapping->fact1;
	double fact2 = mapping->fact2;
	double a = (mapping->L_carac/2.0);
	double eps = (mapping->epsOa)*a;
	double delta = mapping->delta;
	double gamma = (mapping->AoA/180)*M_PI;
	double b = mapping->bj;
	double r0 = mapping->reg*a;

	double t = xi1;
	double theta = xi2;
	double r = r0 + fact2*(exp(fact1*xi1) -1);

	if(x!=NULL){
		*x = ((r*r*r)*cos(theta-gamma) - eps*r*r*cos(delta+gamma) + b*b*r*cos(gamma+theta) - r*r*eps*cos(2*theta-gamma+delta) \
		 + eps*eps*r*cos(theta-gamma) - b*b*eps*cos(delta-gamma))/(r*r + eps*eps - 2.0*r*eps*cos(theta+delta)) + 0.0; //0.0 = shiftx
	}

	if(y!=NULL){
    	*y = ((r*r*r)*sin(theta-gamma) + eps*r*r*sin(delta+gamma) - b*b*r*sin(gamma+theta) - r*r*eps*sin(2*theta-gamma+delta) \
         + eps*eps*r*sin(theta-gamma) - eps*b*b*sin(delta-gamma))/(r*r + eps*eps - 2.0*r*eps*cos(theta+delta)) + 0.0; //0.0 = shifty
	}

	if(h1!=NULL){
		*h1 = sqrt((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+ \
		  (b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)\
		  -(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-\
		  eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)\
		  +(eps*eps)*r*sin(gamma-theta)),2.0)+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*cos(gamma-theta)+\
		  (r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)\
		  /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-\
		  eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-\
		  eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)),2.0));
	}

	if(h2!=NULL){
		*h2 = sqrt(pow(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))\
		  /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)\
		  -eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-\
		  theta))*2.0,2.0)+pow(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))\
		  /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-\
		  eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-\
		  theta))*2.0,2.0));
	}

	if(dh1_dxi1!=NULL){
    	*dh1_dxi1 = fact1*fact2*exp(fact1*t)*1.0/sqrt((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0\
			  +(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0\
			  -eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+\
			  (b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)),2.0)+\
			  (fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-\
			  eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0\
			  /pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-\
		      eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)),2.0))*((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)\
		      *(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)\
			  /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)\
			  -eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)))\
			  *(1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)\
			  *sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0-(eps*sin(delta+gamma)*-2.0+eps*sin(delta-gamma+theta*2.0)\
		      *2.0+r*sin(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps+r*r-eps*r*\
			  cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)\
			  *eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*\
			  ((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*2.0)\
			  *2.0+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*\
			  cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0\
			  /pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*\
			  cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)))*((eps*cos(delta+gamma)*2.0+eps*cos(delta-gamma+theta*2.0)\
			  *2.0-r*cos(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)\
			  -eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))\
			  *2.0+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0\
			  +(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*2.0-pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/\
			  pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*\
			  cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0)*2.0)*(-1.0/2.0);
	}

	if(dh1_dxi2!=NULL){
		*dh1_dxi2 = ((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)\
				 *2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*\
				 cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-\
				 (b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)))*(-((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*sin(gamma+theta)+\
				 eps*r*sin(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)\
				 *2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))+\
				 eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)\
				 -eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
				 cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*\
				 cos(delta-gamma+theta*2.0)*2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)\
				 *cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*\
				 cos(gamma-theta))*4.0)*2.0+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)\
				 +eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/\
				 pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)\
				 +(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)))*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-eps*r*\
				 cos(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)\
			     *((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))+eps*sin(delta+theta)*1.0/\
			     pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)\
				 +(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*\
				 sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*2.0-eps*r*sin(delta+theta)\
			     *(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*\
				 sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*4.0)*2.0)*1.0/sqrt((fact1*fact1)*(fact2*fact2)\
				 *exp(fact1*t*2.0)*pow(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*\
			     sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*\
				 ((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*\
				 r*sin(gamma-theta)),2.0)+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-\
				 eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+\
				 r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-\
				 (b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)),2.0))*(-1.0/2.0);
	}

	if(dh2_dxi1!=NULL){
		*dh2_dxi1 = fact1*fact2*exp(fact1*t)*((((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))\
			  /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-\
			  eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0)\
			  *(-((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)\
		      *2.0)+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+\
			  eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)\
			  *cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*\
			  cos(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+\
			  (b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0\
			  /pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*\
			  cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*4.0)*2.0-(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-\
			  eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+\
			  r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)\
			  +(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-\
			  eps*r*cos(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*\
			  cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))+\
			  eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+\
			  eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
		      cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*\
			  sin(delta+gamma)*2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*\
			  sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*\
			  sin(gamma-theta))*4.0)*2.0)*1.0/sqrt(pow(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*\
			  sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*\
			  cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*\
			  cos(gamma-theta))*2.0,2.0)+pow(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))\
			  /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-\
			  eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0,2.0))*(-1.0/2.0);
	}

	if(dh2_dxi2!=NULL){
		*dh2_dxi2 = 1.0/sqrt(pow(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))/\
				  (eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*\
				  cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0,2.0)+\
				  pow(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/(eps*eps+r*r-eps*r*\
				  cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)\
				  +(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0,2.0))*((((r*r*r)*\
				  cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)\
				  +eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)\
				  +eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0)*(((r*r*r)*sin(gamma-theta)+(b*b)*r*sin(gamma+theta)\
				  +eps*(r*r)*sin(delta-gamma+theta*2.0)*4.0+(eps*eps)*r*sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r\
				  -eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)\
			      *eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)\
			      -(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))*4.0-(eps*eps)*(r*r)*pow(sin(delta+theta),2.0)*1.0\
				  /pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma\
				  +theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*8.0)*2.0-(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma\
			      +theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)\
			      *((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*\
			      r*cos(gamma-theta))*2.0)*(((r*r*r)*cos(gamma-theta)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*4.0+(eps*eps)*r*cos(gamma-theta))/\
				  (eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)\
			      *cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0-(eps*eps)*\
			      (r*r)*pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*\
				  cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*8.0+eps*r*sin(delta+theta)*1.0/\
		          pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*\
				  sin(gamma-theta))*4.0)*2.0)*(1.0/2.0);
	}

	if(d2h1_dxi1dxi2!=NULL){
		*d2h1_dxi1dxi2 = fact1*fact2*exp(fact1*t)*1.0/sqrt(pow(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0\
		                 +(eps*eps)*r*sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+\
						 theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0,2.0)+pow(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)\
						 -eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)\
						 *1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+\
						 eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0,2.0))*((((r*r*r)*cos(gamma-theta)+\
						 (b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*4.0+(eps*eps)*r*cos(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)\
						 +eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*\
						 cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0-(eps*eps)*(r*r)*\
						 pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+\
						 (b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*8.0+eps*r*\
						 sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*\
						 sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))*4.0)*(-((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*\
						 sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/\
						 pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0\
						 +(eps*eps)*r*sin(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*\
						 cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0+\
						 eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*\
						 cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)\
						 *1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*\
						 cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*4.0)*2.0+(((r*r*r)*sin(gamma-theta)+(b*b)*r*sin(gamma+theta)\
						 +eps*(r*r)*sin(delta-gamma+theta*2.0)*4.0+(eps*eps)*r*sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*cos(delta+theta)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*\
						 sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))\
						 *4.0-(eps*eps)*(r*r)*pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*\
						 sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*8.0)\
						 *(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-eps*r*cos(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-(b*b)*\
						 r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)\
						 +(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*\
						 ((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)\
						 *2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps\
						 *(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))\
						 *4.0)*2.0-(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/\
						 (eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)\
						 -eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))\
						 *2.0)*(-((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*8.0)/(eps*eps+r*r-eps\
						 *r*cos(delta+theta)*2.0)+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)+(b*b)*\
						 r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*4.0+(eps*eps)*r*sin(gamma-theta))-eps*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+\
						 (b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)\
						 *cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))*4.0+(eps*eps)*r*\
						 pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*\
						 sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*1.6E1+eps*r*sin(delta+theta)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-eps*r*\
						 cos(delta-gamma+theta*2.0)*4.0)*4.0-eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)\
						 +(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*2.0+(eps*eps)*(r*r)*\
						 pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*\
						 sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*8.0+eps*r*cos(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*\
						 sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*4.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0\
						 +(eps*eps)*r*cos(gamma-theta))*8.0-(eps*eps)*(r*r)*pow(sin(delta+theta),2.0)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,4.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)\
						 +(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.4E1)*2.0-(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*\
						 sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/\
						 pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*\
						 cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)\
						 *3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta-gamma+theta*2.0)*8.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*\
						 4.0+(eps*eps)*r*cos(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*\
						 sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))*4.0+eps*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-\
						 (b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0-(eps*eps)*r*pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*\
						 2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta\
						 -gamma)+(eps*eps)*r*cos(gamma-theta))*1.6E1+eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma\
						 -theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*2.0-(eps*eps)*\
						 (r*r)*pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+\
						 (b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*8.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r\
						 *cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)\
						 *4.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-\
						 (b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))*8.0-eps*r*cos(delta+theta)*(r*2.0-eps*cos(delta+\
						 theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-\
						 eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*4.0+(eps*eps)*(r*r)*pow(sin(delta+theta),2.0)*\
						 (r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,4.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+\
						 (b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.4E1)*2.0)*(1.0/2.0)+\
						 fact1*fact2*exp(fact1*t)*((((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma\
						 -theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)\
						 -eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))\
						 *2.0)*(-((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0)+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*\
						 sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)\
						 -(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*\
						 ((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)\
						 *2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-\
						 eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))\
						 *4.0)*2.0-(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/\
						 (eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-\
						 eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))\
						 *2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-eps*r*cos(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*\
						 cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)\
						 +(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*\
						 ((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)\
						 *2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-\
						 eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))\
						 *4.0)*2.0)*1.0/pow(pow(((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))\
						 /(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-\
						 eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))\
						 *2.0,2.0)+pow(((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))/\
						 (eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)\
						 -eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*\
						 sin(gamma-theta))*2.0,2.0),3.0/2.0)*((((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+\
						 (eps*eps)*r*cos(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*\
						 2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta\
						 -gamma)+(eps*eps)*r*sin(gamma-theta))*2.0)*(((r*r*r)*sin(gamma-theta)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*4.0+\
						 (eps*eps)*r*sin(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*\
						 2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*\
						 sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*\
						 cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))*4.0-(eps*eps)*(r*r)*\
						 pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+\
						 (b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*8.0)*2.0-\
						 (((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))/\
						 (eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*\
						 cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)\
						 +(eps*eps)*r*cos(gamma-theta))*2.0)*(((r*r*r)*cos(gamma-theta)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*4.0+(eps*eps)\
						 *r*cos(gamma-theta))/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+eps*r*cos(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)\
						 *((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*\
						 cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0-(eps*eps)*(r*r)*pow(sin(delta+theta),2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)\
						 *2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*\
						 cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*8.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*\
						 sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))*4.0)*2.0)*(1.0/4.0);
	}

	if(d2h2_dxi1dxi2!=NULL){
		*d2h2_dxi1dxi2 = fact1*fact2*exp(fact1*t)*1.0/sqrt((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)\
						 *3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)\
						 -(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)\
						 +(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)),2.0)+(fact1*fact1)\
						 *(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*\
						 2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-\
						 eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-\
						 gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)),2.0))*((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(1.0/\
						 pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)\
						 *sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0-(eps*sin(delta+gamma)*-2.0+eps*sin(delta-gamma\
						 +theta*2.0)*2.0+r*sin(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps\
						 +r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-\
						 gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*\
						 r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)\
						 *2.0-eps*r*sin(delta+gamma)*2.0)*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-eps*r*cos(delta-gamma\
						 +theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)\
						 *((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))+eps*sin(delta+theta)\
						 *1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)\
						 *sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-\
						 eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+\
						 theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+\
						 theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*\
						 sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*4.0)*2.0+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*sin(gamma-theta)+(r*r)*\
						 sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+\
						 theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*\
						 sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)))*\
						 (1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)\
						 *2.0+(eps*eps)*r*cos(gamma-theta))*2.0+(eps*cos(delta-gamma+theta*2.0)*4.0-r*cos(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-\
						 pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-\
						 eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*cos(gamma-theta))*2.0+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-eps*r*cos(delta-gamma+theta*2.0)*4.0)*2.0\
						 -eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*\
						 sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*4.0-eps*r*sin(delta+theta)*(eps*sin(delta+gamma)*-2.0+\
						 eps*sin(delta-gamma+theta*2.0)*2.0+r*sin(gamma-theta)*6.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*2.0+eps*r*sin(delta+theta)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)\
						 *sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*8.0+eps*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*\
						 2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+\
						 eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*8.0+eps*r*sin(delta+theta)*(r*2.0-eps*\
						 cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*\
						 sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*8.0-eps*r*sin(delta+theta)*pow(r*2.0-eps*cos(delta+theta)*\
						 2.0,2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,4.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+\
						 eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*1.2E1)*2.0+(fact1*fact1)*(fact2*fact2)*\
						 exp(fact1*t*2.0)*((eps*cos(delta+gamma)*2.0+eps*cos(delta-gamma+theta*2.0)*2.0-r*cos(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*\
						 2.0)+1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-\
						 eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/\
						 pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*\
						 cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*2.0-pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+\
						 theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0)*(-((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0-\
						 (b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+(r*2.0-eps*cos(delta+theta)*\
						 2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-\
						 gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*\
						 r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-\
						 gamma)+(eps*eps)*r*cos(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*\
						 cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*\
						 2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*\
						 cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-\
						 gamma)+(eps*eps)*r*cos(gamma-theta))*4.0)*2.0+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)\
						 *cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-\
						 eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*\
						 cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-\
						 gamma)+(eps*eps)*r*cos(gamma-theta)))*((eps*sin(delta-gamma+theta*2.0)*4.0+r*sin(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+\
						 theta)*2.0)-1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*\
						 sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))*2.0+pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps+r*r-eps*\
						 r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*\
						 r*sin(gamma-theta))*2.0-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-\
						 theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)*2.0-eps*sin(delta+theta)*1.0/\
						 pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-\
						 eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*4.0+eps*r*sin(delta+theta)*(eps*cos(delta+gamma)*2.0+eps*\
						 cos(delta-gamma+theta*2.0)*2.0-r*cos(gamma-theta)*6.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*2.0+eps*r*\
						 sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+\
						 (b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*8.0+\
						 eps*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-\
						 theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+\
						 (eps*eps)*r*cos(gamma-theta))*8.0-eps*r*sin(delta+theta)*pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,4.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-\
						 gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*1.2E1+eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+\
						 theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)\
						 *cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*8.0)*2.0)*(1.0/2.0)-fact1*fact2*exp(fact1*t)\
						 *((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-\
						 eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+\
						 theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*\
						 cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)))*(-((eps*eps)*sin(gamma\
						 -theta)+(r*r)*sin(gamma-theta)*3.0-(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)\
						 *2.0)+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-(b*b)*r*sin(gamma\
						 +theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)*2.0+(eps*eps)*r*sin(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*\
						 cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+\
						 theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+\
						 theta)*2.0,2.0)*((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*\
						 cos(delta-gamma+theta*2.0)*2.0)*2.0-eps*r*sin(delta+theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+\
						 theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)\
						 -(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta))*4.0)*2.0+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*\
						 sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*\
						 2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)\
						 *((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*\
						 sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)))*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0-(b*b)*cos(gamma+theta)-\
						 eps*r*cos(delta-gamma+theta*2.0)*4.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+\
						 r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)*2.0+\
						 (eps*eps)*r*cos(gamma-theta))+eps*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-\
						 eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*\
						 sin(gamma-theta))*2.0+eps*r*sin(delta+theta)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*\
						 sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*2.0-eps*r*sin(delta+\
						 theta)*(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*\
						 sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-\
						 theta))*4.0)*2.0)*1.0/pow((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+\
						 (b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-\
						 (r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+\
						 gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)),2.0)+\
						 (fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*pow(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-\
						 eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+\
						 theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*\
						 cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)),2.0),3.0/2.0)*\
						 ((fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+\
						 eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta\
						 +theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*\
						 sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)))*(1.0/pow(eps*eps+\
						 r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*\
						 sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0-(eps*sin(delta+gamma)*-2.0+eps*sin(delta\
						 -gamma+theta*2.0)*2.0+r*sin(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*\
						 1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)\
						 +eps*(r*r)*sin(delta-gamma+theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta))*2.0+(r*2.0-eps*cos(delta+theta)*\
						 2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+\
						 theta)+eps*r*sin(delta-gamma+theta*2.0)*2.0-eps*r*sin(delta+gamma)*2.0)*2.0)*2.0+(fact1*fact1)*(fact2*fact2)*exp(fact1*t*2.0)*\
						 (((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+\
						 theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+\
						 theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*\
						 2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)))*((eps*cos(delta+gamma)*2.0+eps*cos(delta-gamma+theta*2.0)*2.0-r*\
						 cos(gamma-theta)*6.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)+1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*\
						 cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)\
						 +(eps*eps)*r*cos(gamma-theta))*2.0+(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,2.0)*((eps*eps)*\
						 cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-eps*r*cos(delta-gamma+theta*2.0)*2.0)*\
						 2.0-pow(r*2.0-eps*cos(delta+theta)*2.0,2.0)*1.0/pow(eps*eps+r*r-eps*r*cos(delta+theta)*2.0,3.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*\
						 cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-\
						 theta))*2.0)*2.0)*(1.0/4.0);
	}

	double dxdxi1 = fact1*fact2*exp(fact1*t)*(((eps*eps)*cos(gamma-theta)+(r*r)*cos(gamma-theta)*3.0+(b*b)*cos(gamma+theta)-eps*r*cos(delta+gamma)*2.0-\
			     eps*r*cos(delta-gamma+theta*2.0)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-\
				 eps*r*cos(delta+theta)*2.0,2.0)*((r*r*r)*cos(gamma-theta)-eps*(r*r)*cos(delta+gamma)+(b*b)*r*cos(gamma+theta)-eps*(r*r)*cos(delta-\
				 gamma+theta*2.0)-(b*b)*eps*cos(delta-gamma)+(eps*eps)*r*cos(gamma-theta)));

    double dydxi1 = -fact1*fact2*exp(fact1*t)*(((eps*eps)*sin(gamma-theta)+(r*r)*sin(gamma-theta)*3.0+(b*b)*sin(gamma+theta)+eps*r*sin(delta-gamma+theta*2.0)\
				  *2.0-eps*r*sin(delta+gamma)*2.0)/(eps*eps+r*r-eps*r*cos(delta+theta)*2.0)-(r*2.0-eps*cos(delta+theta)*2.0)*1.0/pow(eps*eps+r*r-eps*r*\
				  cos(delta+theta)*2.0,2.0)*((r*r*r)*sin(gamma-theta)-eps*(r*r)*sin(delta+gamma)+(b*b)*r*sin(gamma+theta)+eps*(r*r)*sin(delta-gamma+\
				  theta*2.0)+(b*b)*eps*sin(delta-gamma)+(eps*eps)*r*sin(gamma-theta)));

	if(theta_mesh!=NULL){
    	*theta_mesh = atan2(dydxi1,dxdxi1);
	}

}

void airfoil_init_velocity(AirfoilMapping* mapping, double U_inf, double xi1, double xi2, double * u_n, double *u_theta){

	std::complex<double> I(0,1);

	double a = (mapping->L_carac/2.0);
    double r0 = mapping->reg*a;
	double fact1 = mapping->fact1;
	double fact2 = mapping->fact2;
	double eps = (mapping->epsOa)*a;
	double delta = mapping->delta;
	double gamma = (mapping->AoA/180)*M_PI;
	double b = mapping->bj;
    double Beta_tr = std::real(-std::log((eps*std::exp(-I*delta)+b)/a)/I);
    double Gamma = 0.0;

	//Value of theta on the mesh :
	double theta_mesh;
	airfoil_metrics(mapping,xi1,xi2, NULL, NULL, \
							  NULL, NULL, NULL, NULL,\
							  NULL, NULL, \
							  NULL, NULL, \
							  &theta_mesh);

    double r;
	std::complex<double> z, vel;

	//u_n component :
	if(u_n!=NULL){
		r = r0 + fact2*exp(xi1*fact1) - fact2;
		z = r*std::exp(I * xi2);

		/*be careful vel = u_x - i*u_y*/
		vel = ((U_inf*(std::exp(-I*gamma)-((a*a)*std::exp(I*gamma))/(std::pow(z,2.0))) + Gamma/(2.0*M_PI*I*z) )\
		/(1.0 - b*b/(std::pow(z-eps*std::exp(-I*delta),2.0))))*std::exp(I*gamma);

		*u_n = std::real(vel)*cos(theta_mesh) - std::imag(vel)*sin(theta_mesh);
	}


	//u_theta component :
	if(u_theta!=NULL){
		r = r0+ fact2*exp( xi1*fact1 ) - fact2;
		z = r*std::exp(I * xi2);

        /*be careful vel = u_x - i*u_y*/
        vel = ((U_inf*(std::exp(-I*gamma)-((a*a)*std::exp(I*gamma))/(std::pow(z,2.0))) + Gamma/(2.0*M_PI*I*z) )\
                  /(1.0 - b*b/(std::pow(z-eps*std::exp(-I*delta),2.0))))*std::exp(I*gamma);

		*u_theta = -std::real(vel)*sin(theta_mesh) - std::imag(vel)*cos(theta_mesh);
	}
}

void airfoil_init_pressure(AirfoilMapping* mapping, double U_inf, double xi1, double xi2, double * P){

	std::complex<double> I(0,1);

	double a = (mapping->L_carac/2.0);
    double r0 = mapping->reg*a;
	double fact1 = mapping->fact1;
	double fact2 = mapping->fact2;
	double eps = (mapping->epsOa)*a;
	double delta = mapping->delta;
	double gamma = (mapping->AoA/180)*M_PI;
	double b = mapping->bj;
 	double Beta_tr = std::real(-std::log((eps*std::exp(-I*delta)+b)/a)/I);
 	double Gamma = 0.0;

 	double r;
	std::complex<double> z, vel;

 	r = r0 + fact2*exp(xi1*fact1) - fact2;

	z = r*std::exp(I * xi2);

 	/*be careful vel = u_x - i*u_y*/
 	vel = ((U_inf*(std::exp(-I*gamma)-((a*a)*std::exp(I*gamma))/(std::pow(z,2.0))) + Gamma/(2.0*M_PI*I*z) )\
            /(1.0 - b*b/(std::pow(z-eps*std::exp(-I*delta),2.0))))*std::exp(I*gamma);

	*P = -(std::real(vel)*std::real(vel) + std::imag(vel)*std::imag(vel))/2.0;
}
