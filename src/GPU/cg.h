#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <type_traits>

cublasStatus_t cublas_axpy(cublasHandle_t handle, int n, const double *alpha, const double *x, int incx, double * y, int incy )
{ return cublasDaxpy(handle, n, alpha, x, incx, y, incy); }

cublasStatus_t cublas_axpy(cublasHandle_t handle, int n, const float *alpha, const float *x, int incx, float * y, int incy )
{ return cublasSaxpy(handle, n, alpha, x, incx, y, incy); }

cublasStatus_t cublas_dot(cublasHandle_t handle, int n, const double *x, int incx, const double * y, int incy, double *result )
{ return cublasDdot(handle,n,x,incx,y,incy,result); }

cublasStatus_t cublas_dot(cublasHandle_t handle, int n, const float *x, int incx, const float * y, int incy, float *result )
{ return cublasSdot(handle,n,x,incx,y,incy,result); }

cublasStatus_t cublas_scal(cublasHandle_t handle, int n, const float *x, float * y, int incy )
{ return cublasSscal(handle,n,x,y,incy); }

cublasStatus_t cublas_scal(cublasHandle_t handle, int n, const double *x, double * y, int incy )
{ return cublasDscal(handle,n,x,y,incy); }

cublasStatus_t cublas_copy(cublasHandle_t handle, int n, const double *x, int incx, double * y, int incy)
{ return cublasDcopy(handle,n,x,incx,y,incy); }

cublasStatus_t cublas_copy(cublasHandle_t handle, int n, const float *x, int incx, float * y, int incy)
{ return cublasScopy(handle,n,x,incx,y,incy); }


namespace GPU
{
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

template <typename T>
T cg(int *I, int *J, T *val, T *x, T *rhs,const int N,const int nz,const T tol, const int max_iter, int &nb_iter)
{
cudaDataType_t size_float;

if (std::is_same<T,float>::value)
	{ size_float = CUDA_R_32F; }
else if (std::is_same<T,double>::value)
	{ size_float = CUDA_R_64F; }
else exit(1);

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

  /* Wrap raw data into cuSPARSE generic API objects */
cusparseSpMatDescr_t matA;
cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val,CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO, size_float);

//cudaMemcpy(d_x,x,N*sizeof(T),cudaMemcpyHostToDevice);// bug !!

cusparseDnVecDescr_t vecx;
cusparseCreateDnVec(&vecx, N, d_x, size_float);

cusparseDnVecDescr_t vecp;
cusparseCreateDnVec(&vecp, N, d_p, size_float);

cusparseDnVecDescr_t vecAx;
cusparseCreateDnVec(&vecAx, N, d_Ax, size_float);

cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_val, val, nz * sizeof(T), cudaMemcpyHostToDevice);
cudaMemcpy(d_x,x,N*sizeof(T), cudaMemcpyHostToDevice);
cudaMemcpy(d_r,rhs,N*sizeof(T), cudaMemcpyHostToDevice);

T alpha = 1.0;
T alpham1 = -1.0;
T beta = 0.0;
r0 = 0.0;

/* Allocate workspace for cuSPARSE */
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

cudaMemcpy(x, d_x, N * sizeof(T), cudaMemcpyDeviceToHost);

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

} // end namespace GPU


