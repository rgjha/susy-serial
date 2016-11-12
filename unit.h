#include "utilities.h"

void unit(const Gauge_Field &U, Gauge_Field &U2);
extern "C" void zgeev_( char*, char*, int*, double at[], int *, double b[], double dummy[], 
            int *, double dummy2[], int*, double work[], int *, double work2[], int *);
