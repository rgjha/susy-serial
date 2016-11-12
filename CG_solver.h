#include "utilities.h"
#include "matrix.h"

#ifndef CG_RESIDUAL_H
#define CG_RESIDUAL_H
const double CG_RESIDUAL = 0.0000000001;
#endif


void CG_solver(const Gauge_Field &U, 
	       const Twist_Fermion &rhs, Twist_Fermion &sol);
