double *vector(int,int);
double **matrix(int,int,int,int);
void free_vector();
void free_matrix();

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double wr[], double wi[]);
