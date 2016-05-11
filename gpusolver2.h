#include "utilities.h"

// Calls gpu routines for solving the linear system
// the solution is left in the pointer soln
void gpusolver2(Complex *m, int *col, int *row,
    Complex *bn, Complex *sol);
