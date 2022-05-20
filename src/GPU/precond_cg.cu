#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <type_traits>

#include "../alg_utils.h"
//#include "gpu_utils.h"

/**
conjuguate gradient with diagonal preconditionner

*/

template <class T>
T _precond_cg(alg::CSR_mat<T> const& A, T *x, T *rhs, const T tol, const int max_iter, int &nb_iter)
{
int k(0);
T r(0);
cudaDataType_t size_float;

if (std::is_same<T,float>::value)
	{ size_float = CUDA_R_32F; }
else if (std::is_same<T,double>::value)
	{ size_float = CUDA_R_64F; }
else exit(1);

const int nz = A.I[A.N];//nb of coefficients of the sparse matrix
const int N = A.N;

cublasHandle_t cublasHandle = 0;
cublasCreate(&cublasHandle);

cusparseHandle_t cusparseHandle = 0;
cusparseCreate(&cusparseHandle);



cusparseDestroy(cusparseHandle);
cublasDestroy(cublasHandle);

nb_iter = k;
return sqrt(r);
}

double precond_cg(alg::CSR_mat<double> const& A, double *x, double *rhs, const double tol, const int max_iter, int &nb_iter)
{ return _precond_cg<double>(A, x, rhs, tol, max_iter, nb_iter); }

float precond_cg(alg::CSR_mat<float> const& A, float *x, float *rhs, const float tol, const int max_iter, int &nb_iter)
{ return _precond_cg<float>(A, x, rhs, tol, max_iter, nb_iter); }
