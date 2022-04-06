#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include<iostream>

/*
nz = nb de coefficients, stockés dans val, I,J indirections sur les indices de la matrices creuses x résultat ; rhs, second membre du système, de dimension N
*/

void cg(int *I, int *J, float *val, float *x, float *rhs, int N, int nz,const float tol, int max_iter)
{
cudaDataType_t size_float = CUDA_R_32F;

int *d_col, *d_row;
int k;
float *d_val, *d_x, *d_r, *d_p, *d_Ax;
float r0,r1,a,na,b,dot;

cublasHandle_t cublasHandle = 0;
cublasCreate(&cublasHandle);

cusparseHandle_t cusparseHandle = 0;
cusparseCreate(&cusparseHandle);


cudaMalloc((void **)&d_col, nz * sizeof(int));
cudaMalloc((void **)&d_row, (N + 1) * sizeof(int));
cudaMalloc((void **)&d_val, nz * sizeof(float));

cudaMalloc((void **)&d_x, N * sizeof(float));
cudaMalloc((void **)&d_r, N * sizeof(float));
cudaMalloc((void **)&d_p, N * sizeof(float));
cudaMalloc((void **)&d_Ax, N * sizeof(float));

  /* Wrap raw data into cuSPARSE generic API objects */
cusparseSpMatDescr_t matA;
cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val,CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO, size_float);

cudaMemcpy(d_x,x,N*sizeof(float),cudaMemcpyHostToDevice);

cusparseDnVecDescr_t vecx;
cusparseCreateDnVec(&vecx, N, d_x, size_float);

cusparseDnVecDescr_t vecp;
cusparseCreateDnVec(&vecp, N, d_p, size_float);

cusparseDnVecDescr_t vecAx;
cusparseCreateDnVec(&vecAx, N, d_Ax, size_float);

cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(d_x,x,N*sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(d_r,rhs,N*sizeof(float), cudaMemcpyHostToDevice);

float alpha = 1.0;
float alpham1 = -1.0;
float beta = 0.0;
r0 = 0.0;

/* Allocate workspace for cuSPARSE */
size_t bufferSize = 0;
cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx,&beta, vecAx, size_float, CUSPARSE_CSRMV_ALG1, &bufferSize);
void *buffer = NULL;
cudaMalloc(&buffer, bufferSize);

cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, size_float,CUSPARSE_CSRMV_ALG1, buffer);

cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
cublasSdot(cublasHandle, N, d_r , 1, d_r, 1, &r1);

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
	
	cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecp, &beta, vecAx, size_float, CUSPARSE_CSRMV_ALG1, buffer);//CUSPARSE_SPMV_ALG_DEFAULT
	
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

cudaMemcpy(x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);

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

void cg(int *I, int *J, double *val, double *x, double *rhs, int N, int nz,const double tol, int max_iter)
{
cudaDataType_t size_float = CUDA_R_64F;

int *d_col, *d_row;
int k;
double *d_val, *d_x, *d_r, *d_p, *d_Ax;
double r0,r1,a,na,b,dot;

cublasHandle_t cublasHandle = 0;
cublasCreate(&cublasHandle);

cusparseHandle_t cusparseHandle = 0;
cusparseCreate(&cusparseHandle);


cudaMalloc((void **)&d_col, nz * sizeof(int));
cudaMalloc((void **)&d_row, (N + 1) * sizeof(int));
cudaMalloc((void **)&d_val, nz * sizeof(double));

cudaMalloc((void **)&d_x, N * sizeof(double));
cudaMalloc((void **)&d_r, N * sizeof(double));
cudaMalloc((void **)&d_p, N * sizeof(double));
cudaMalloc((void **)&d_Ax, N * sizeof(double));

  /* Wrap raw data into cuSPARSE generic API objects */
cusparseSpMatDescr_t matA;
cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val,CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO, size_float);

cudaMemcpy(d_x,x,N*sizeof(double),cudaMemcpyHostToDevice);

cusparseDnVecDescr_t vecx;
cusparseCreateDnVec(&vecx, N, d_x, size_float);

cusparseDnVecDescr_t vecp;
cusparseCreateDnVec(&vecp, N, d_p, size_float);

cusparseDnVecDescr_t vecAx;
cusparseCreateDnVec(&vecAx, N, d_Ax, size_float);

cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_val, val, nz * sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(d_x,x,N*sizeof(double), cudaMemcpyHostToDevice);
cudaMemcpy(d_r,rhs,N*sizeof(double), cudaMemcpyHostToDevice);

double alpha = 1.0;
double alpham1 = -1.0;
double beta = 0.0;
r0 = 0.0;

/* Allocate workspace for cuSPARSE */
size_t bufferSize = 0;
cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx,&beta, vecAx, size_float, CUSPARSE_CSRMV_ALG1, &bufferSize);
void *buffer = NULL;
cudaMalloc(&buffer, bufferSize);

cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, size_float,CUSPARSE_CSRMV_ALG1, buffer);

cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
cublasDdot(cublasHandle, N, d_r , 1, d_r, 1, &r1);


cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
cublasDdot(cublasHandle, N, d_r , 1, d_r, 1, &r1);

k=1;

while(r1 > tol*tol && k <= max_iter)
	{
	if (k>1)
		{
		b = r1/r0;
		cublasDscal(cublasHandle, N, &b, d_p, 1);
		cublasDaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
		}
	else
		{ cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1); }
	
	cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecp, &beta, vecAx, size_float, CUSPARSE_CSRMV_ALG1, buffer);//CUSPARSE_SPMV_ALG_DEFAULT
	
	cublasDdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
	a = r1/dot;
	
	cublasDaxpy(cublasHandle,N,&a,d_p,1,d_x,1);
	na = -a;
	cublasDaxpy(cublasHandle,N,&na,d_Ax,1,d_r,1);
	
	r0 = r1;
	cublasDdot(cublasHandle,N,d_r,1,d_r,1,&r1);
	cudaDeviceSynchronize();
	std::cout << "iteration = " << k << " ; residual = " << sqrt(r1) << std::endl; 
	k++;
	}

cudaMemcpy(x, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);

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

