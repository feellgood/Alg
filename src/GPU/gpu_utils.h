#ifndef GPU_UTILS_H
#define GPU_UTILS_H

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

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

#endif
