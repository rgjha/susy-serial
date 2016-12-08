#include "utilities.h"
int lat_pack(const Lattice_Vector &x);
void lat_unpack(const int n, Lattice_Vector &x);
int pack(const int flavor,const int type, const Lattice_Vector &x, const int color);
void unpack(const int i, int &flavor,int &type, Lattice_Vector &x, int &color);

#ifdef FULLMATRIX
void full_fermion_op(const Gauge_Field &U,Complex M[LEN][LEN]);
Complex Pfaffian(Complex M[LEN][LEN]);
void eigenvalues(Complex M[LEN][LEN]);
#endif
double cnorm(const Complex &c);
void test(const Gauge_Field &);

extern "C" void zgeev_( char*, char*, int*, double at[], int *, double b[], double dummy[],
            int *, double dummy2[], int*, double work[], int *, double work2[], int *);

