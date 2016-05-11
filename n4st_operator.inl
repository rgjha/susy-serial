#include <cusp/linear_operator.h>

// to use this: in gpusolver.cu
// #include this file
//
// then once you have formed csr matrices for M, M^\\dagger, say they're called
// M, M_H, insert the line:
// dirac_operator AA(M,M_H);
// then in the cg_m call use AA as the matrix

// typedefs for short hand
typedef cusp::complex<double> ScalarType;
typedef cusp::device_memory MemorySpace;
typedef cusp::csr_matrix<int,ScalarType,MemorySpace> Matrix;

// the clase we want to create
// the constructor takes two csr matrices
class dirac_operator :
  public cusp::linear_operator<ScalarType,MemorySpace>
{
    public:
    // all examples have this... not exactly sure why... i think it is so that
    // if foo is an instance of this class then foo.num_rows is correct.
    typedef cusp::linear_operator<ScalarType,MemorySpace> super;

    // matrix
    Matrix M;
    // its adjoint
    Matrix M_H;

    size_t num_rows;
    size_t num_cols;
    
    // constructor
    dirac_operator(Matrix _M, Matrix _M_H) :
                  super(M.num_cols,M.num_rows), M(_M), M_H(_M_H), num_rows(_M_H.num_rows), num_cols(_M.num_cols) {}

    // linear operator y = A*x... needs to be this way so that cg or
    // cg_m can call cusp::multiply on it
    template <typename VectorType1,
              typename VectorType2>
    void operator()(const VectorType1& x, VectorType2& y) const
    {
      // temporary variable
      VectorType1 xx(x.size());
      // do the first multiplication and write to xx
      cusp::multiply(M,x,xx);
      // then do the second multiplication and write to y
      cusp::multiply(M_H,xx,y);
      cusp::blas::axpy(x,y,MASS*MASS);
    }
};
