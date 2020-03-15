//
//  thomas.h
//  MECA2660_hw 2016
//
//  Created by Thomas Gillis on 16/02/16.
//  Copyright © 2016 Thomas Gillis. All rights reserved.
//

#ifndef _THOMAS_H_
#define _THOMAS_H_

#include <stdlib.h>

void solve_Ac_thomas(const int n, const double a, const double b, const double c, double* x,double* q);
void solve_period_3diag(const int n, const double a, const double b, const double c, double* x,double* q);

#endif
