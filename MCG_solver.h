#include "utilities.h"
#include "matrix.h"


#ifndef CG_RESIDUAL_H
#define CG_RESIDUAL_H
const double CG_RESIDUAL = 0.0000000001;
#endif


void MCG_solver(const Gauge_Field &, const Twist_Fermion &rhs, double shift[], 
Twist_Fermion sol[DEGREE], 
Twist_Fermion psol[DEGREE]);
