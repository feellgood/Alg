#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <type_traits>

#include "../alg_utils.h"
#include "gpu_utils.h"


/**
I,J indices of the non zero coefficients of the sparse matrix
val : values of the non zero coefficients of the sparse matrix 
x : result of dimension N
rhs : right hand size
N : dimension of the linear system to solve
nz : nb of coefficients of the sparse matrix
tol : tolerance
max_iter : nb iteration maximum
nb_iter : number of iteration done
returns residue
*/

template <class T>
T _cg(alg::CSR_mat<T> const& A, T *x, T *rhs, const T tol, const int max_iter, int &nb_iter)
{
cudaDataType_t size_float;

if (std::is_same<T,float>::value)
	{ size_float = CUDA_R_32F; }
else if (std::is_same<T,double>::value)
	{ size_float = CUDA_R_64F; }
else exit(1);

const int nz = A.I[A.N];
const int N = A.N;

int *d_col, *d_row;
int k;
T *d_val, *d_x, *d_r, *d_p, *d_Ax;
T r0,r1,a,na,b,dot;

cublasHandle_t cublasHandle = 0;
cublasCreate(&cublasHandle);

cusparseHandle_t cusparseHandle = 0;
cusparseCreate(&cusparseHandle);

cudaMalloc((void **)&d_col, nz * sizeof(int));
cudaMalloc((void **)&d_row, (N + 1) * sizeof(int));
cudaMalloc((void **)&d_val, nz * sizeof(T));

cudaMalloc((void **)&d_x, N * sizeof(T));
cudaMalloc((void **)&d_r, N * sizeof(T));
cudaMalloc((void **)&d_p, N * sizeof(T));
cudaMalloc((void **)&d_Ax, N * sizeof(T));

  // Wrap raw data into cuSPARSE generic API objects
cusparseSpMatDescr_t matA;
cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val,CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO, size_float);

cusparseDnVecDescr_t vecx;
cusparseCreateDnVec(&vecx, N, d_x, size_float);

cusparseDnVecDescr_t vecp;
cusparseCreateDnVec(&vecp, N, d_p, size_float);

cusparseDnVecDescr_t vecAx;
cusparseCreateDnVec(&vecAx, N, d_Ax, size_float);

cudaMemcpy(d_col, A.J, nz * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_row, A.I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_val, A.val, nz * sizeof(T), cudaMemcpyHostToDevice);
cudaMemcpy(d_x,x,N*sizeof(T), cudaMemcpyHostToDevice);
cudaMemcpy(d_r,rhs,N*sizeof(T), cudaMemcpyHostToDevice);

T alpha = 1.0;
T alpham1 = -1.0;
T beta = 0.0;
r0 = 0.0;

// Allocate workspace for cuSPARSE 
size_t bufferSize = 0;
cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx,&beta, vecAx, size_float, CUSPARSE_CSRMV_ALG1, &bufferSize);
void *buffer = NULL;
cudaMalloc(&buffer, bufferSize);

cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, size_float,CUSPARSE_CSRMV_ALG1, buffer);

cublas_axpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
cublas_dot(cublasHandle, N, d_r , 1, d_r, 1, &r1);

k=1;

while(r1 > tol*tol && k <= max_iter)
	{
	if (k>1)
		{
		b = r1/r0;
		cublas_scal(cublasHandle, N, &b, d_p, 1);
		cublas_axpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
		}
	else
		{ cublas_copy(cublasHandle, N, d_r, 1, d_p, 1); }
	
	cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecp, &beta, vecAx, size_float, CUSPARSE_CSRMV_ALG1, buffer);//CUSPARSE_SPMV_ALG_DEFAULT
	
	cublas_dot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
	a = r1/dot;
	
	cublas_axpy(cublasHandle,N,&a,d_p,1,d_x,1);
	na = -a;
	cublas_axpy(cublasHandle,N,&na,d_Ax,1,d_r,1);
	
	r0 = r1;
	cublas_dot(cublasHandle,N,d_r,1,d_r,1,&r1);
	cudaDeviceSynchronize();
	k++;
	}

cudaMemcpy(x, d_x, N*sizeof(T), cudaMemcpyDeviceToHost);

cusparseDestroy(cusparseHandle);
cublasDestroy(cublasHandle);

cusparseDestroySpMat(matA);
cusparseDestroyDnVec(vecx);
cusparseDestroyDnVec(vecAx);

cudaFree(d_col);
cudaFree(d_row);
cudaFree(d_val);
cudaFree(d_x);
cudaFree(d_r);
cudaFree(d_p);
cudaFree(d_Ax);

nb_iter = k;
return sqrt(r1);
}

double cg(alg::CSR_mat<double> const& A, double *x, double *rhs, const double tol, const int max_iter, int &nb_iter)
{ return _cg<double>(A, x, rhs, tol, max_iter, nb_iter); }

float cg(alg::CSR_mat<float> const& A, float *x, float *rhs, const float tol, const int max_iter, int &nb_iter)
{ return _cg<float>(A, x, rhs, tol, max_iter, nb_iter); }

