#include "utilities.h"
#include "action.h"

void obs(const Gauge_Field &U, const Twist_Fermion &F,
double &act_s, double &mass, Complex &d,
double &act_F, double eigenvals[SITES][NUMLINK][NCOLOR]);

extern "C" void zgeev_( char*, char*, int*, double at[], int *, double b[], double dummy[], 
            int *, double dummy2[], int*, double work[], int *, double work2[], int *);
