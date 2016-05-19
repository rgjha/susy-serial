#include "gpusolver.h"
#include <cusp/krylov/cg_m.h>
#include <cusp/krylov/cg.h>
#include <cusp/multiply.h>
#include <cusp/print.h>
#include "utilities.h"
#include "n4st_operator.inl"

void gpusolver(Complex *m, int *col, int *row, Complex *bn, double *shift,
    Complex soln[DEGREE][LEN])
{

  // Declare cusp containers on host
<<<<<<< HEAD
  cusp::csr_matrix<int,cusp::complex<double>,cusp::host_memory> A_h(LEN,LEN,LEN*NONZEROES);
  cusp::csr_matrix<int,cusp::complex<double>,cusp::host_memory> A_h_Prime(LEN,LEN,LEN*NONZEROES);
=======
  cusp::csr_matrix<int,cusp::complex<double>,cusp::host_memory> A_h(LEN,LEN,TOTALNONZEROES);
  cusp::csr_matrix<int,cusp::complex<double>,cusp::host_memory> A_h_Prime(LEN,LEN,TOTALNONZEROES);
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
  cusp::array1d<cusp::complex<double>, cusp::host_memory> b_h(LEN,0.0);
  cusp::array1d<double, cusp::host_memory> sigma_h(DEGREE);

  // Copy data into cusp containers
  for (int i=0; i<DEGREE; i++)
    sigma_h[i] = shift[i];
  
  for (int i = 0; i<LEN; i++)
  {
    b_h[i] = cusp::complex<double>(bn[i].real(),bn[i].imag());
    A_h.row_offsets[i] = row[i];
    A_h_Prime.row_offsets[i] = row[i];
    if(row[i] == -1){cout << "-1 detected, exiting in row" << endl; exit(1);}
  }
<<<<<<< HEAD
  for (int i=0; i<LEN*NONZEROES; i++)
=======
  for (int i=0; i<TOTALNONZEROES; i++)
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
  {
    A_h.column_indices[i] = col[i];
    A_h_Prime.column_indices[i] = col[i];
    A_h.values[i] = cusp::complex<double>(m[i].real(),m[i].imag());
    A_h_Prime.values[i] = cusp::complex<double>(-m[i].real(),m[i].imag());
    if(col[i] == -1){cout << "-1 detected, exiting in col" << endl; exit(1);}
  }


<<<<<<< HEAD
  A_h.row_offsets[LEN] = LEN*NONZEROES;
  A_h_Prime.row_offsets[LEN] = LEN*NONZEROES;
=======
  A_h.row_offsets[LEN] = TOTALNONZEROES;
  A_h_Prime.row_offsets[LEN] = TOTALNONZEROES;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
  
  // Copy information to device
  cusp::csr_matrix<int,cusp::complex<double>,cusp::device_memory> A_d = A_h;
  cusp::csr_matrix<int,cusp::complex<double>,cusp::device_memory> A_d_Prime = A_h_Prime;
  cusp::array1d<cusp::complex<double>, cusp::device_memory> b_d = b_h;
  cusp::array1d<double, cusp::device_memory> sigma_d = sigma_h;
  
  // Symmetrize
  dirac_operator AA(A_d, A_d_Prime);

  // Declare the solution on the device
  cusp::array1d<cusp::complex<double>, cusp::device_memory> x_d(LEN*DEGREE, 1);

  // monitor the cg iteration silently
  cusp::default_monitor<cusp::complex<double> > monitor(b_d, 1000, 1e-14);
  
  // Solve
  cusp::krylov::cg_m(AA,x_d,b_d,sigma_d,monitor);
  
  cusp::array1d<cusp::complex<double>, cusp::host_memory> x_h = x_d;
 
  for (int n=0; n<DEGREE; n++)
  {
    for (int i=0; i<LEN; i++)
    {
      soln[n][i] = Complex(x_h[i+LEN*n].real(),x_h[i+LEN*n].imag());
    }
  }
}

