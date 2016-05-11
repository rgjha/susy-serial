#include "utilities.h"
int lat_pack(const Lattice_Vector &x);
void lat_unpack(const int n, Lattice_Vector &x);
int pack(const int flavor,const int type, const Lattice_Vector &x, const int color);
void unpack(const int i, int &flavor,int &type, Lattice_Vector &x, int &color);

void build_sparse_matrix(const Adjoint_Links &V, const Gauge_Field &U, Complex m[LEN], int col[LEN], int row[]);
#ifdef FULLMATRIX
void full_fermion_op(const Gauge_Field &U,Complex M[LEN][LEN]);
Complex Pfaffian(Complex M[LEN][LEN]);
void eigenvalues(Complex M[LEN][LEN]);
#endif
void build_vector(const Twist_Fermion &t, Complex v[LEN]);
Twist_Fermion extract_vector(Complex v[LEN]);
void sparse_mult(Complex m[], int col[], int row[], int, Complex v[], Complex
t[]);
Complex dot(Complex v[], Complex t[]);
double cnorm(const Complex &c);
void test(const Gauge_Field &);

extern "C" void zgeev_( char*, char*, int*, double at[], int *, double b[], double dummy[],
            int *, double dummy2[], int*, double work[], int *, double work2[], int *);

void mysort(Complex d[], int key[], int n);
void col_order(Complex m[], int col[], int row[]);
