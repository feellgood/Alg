#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include<iostream>

/* genTridiag: generate a random tridiagonal symmetric matrix */
template <typename T>
void genTridiag(int *I, int *J, T *val, int N, int nz)
{
  I[0] = 0, J[0] = 0, J[1] = 1;
  val[0] = (T) rand() / RAND_MAX + 10.0f;
  val[1] = (T) rand() / RAND_MAX;
  int start;

  for (int i = 1; i < N; i++) 
  {
    if (i > 1) { I[i] = I[i - 1] + 3; } else { I[1] = 2; }

    start = (i - 1) * 3 + 2;
    J[start] = i - 1;
    J[start + 1] = i;

    if (i < N - 1) { J[start + 2] = i + 1; }

    val[start] = val[start - 1];
    val[start + 1] = (T) rand() / RAND_MAX + 10.0f;

    if (i < N - 1) { val[start + 2] = (T) rand() / RAND_MAX; }
  }

  I[N] = nz;
}

/*
nz = nb de coefficients, stockés dans val, I,J indirections sur les indices de la matrices creuses x résultat ; rhs, second membre du système, de dimension N
*/

template <typename T>
void cg(int *I, int *J, T *val, T *x, T *rhs, int N, int nz, T tol, int max_iter)
{
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
cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val,CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);

cudaMemcpy(d_x,x,N*sizeof(T),cudaMemcpyHostToDevice);

cusparseDnVecDescr_t vecx;
cusparseCreateDnVec(&vecx, N, d_x, CUDA_R_32F);

cusparseDnVecDescr_t vecp;
cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_32F);

cusparseDnVecDescr_t vecAx;
cusparseCreateDnVec(&vecAx, N, d_Ax, CUDA_R_32F);

cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_val, val, nz * sizeof(T), cudaMemcpyHostToDevice);
cudaMemcpy(d_x,x,N*sizeof(T), cudaMemcpyHostToDevice);
cudaMemcpy(d_r,rhs,N*sizeof(T), cudaMemcpyHostToDevice);

T alpha = 1.0f;
T alpham1 = -1.0f;
T beta = 0.0f;
r0 = 0.0f;

/* Allocate workspace for cuSPARSE */
size_t bufferSize = 0;
cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx,&beta, vecAx, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, &bufferSize);
void *buffer = NULL;
cudaMalloc(&buffer, bufferSize);

cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, CUDA_R_32F,CUSPARSE_CSRMV_ALG1, buffer);

cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
cublasSdot(cublasHandle, N, d_r , 1, d_r, 1, &r1);

k=1;

while(r1 > tol*tol && k <= max_iter)
	{
	if (k>1)
		{
		b = r1/r0;
		cublasSscal(cublasHandle, N, &b, d_p, 1);
		cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
		}
	else
		{ cublasScopy(cublasHandle, N, d_r, 1, d_p, 1); }
	
	cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecp, &beta, vecAx, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, buffer);//CUSPARSE_SPMV_ALG_DEFAULT
	
	cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
	a = r1/dot;
	
	cublasSaxpy(cublasHandle,N,&a,d_p,1,d_x,1);
	na = -a;
	cublasSaxpy(cublasHandle,N,&na,d_Ax,1,d_r,1);
	
	r0 = r1;
	cublasSdot(cublasHandle,N,d_r,1,d_r,1,&r1);
	cudaDeviceSynchronize();
	std::cout << "iteration = " << k << " ; residual = " << sqrt(r1) << std::endl; 
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
}
