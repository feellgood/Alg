#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<iostream>

/** template  Kernel to add in parallel on GPU : c = a+b */
template<typename T>
__global__ void addKernel(const T *a, const T *b, T *c, int size)
{
int i = blockIdx.x *blockDim.x + threadIdx.x;
if (i<size)
	{ c[i] = a[i] + b[i]; }
}

template<typename T> 
void cudaAdd(const T *a, const T *b, T *c, int size)
{
T * dev_a = nullptr;
T * dev_b = nullptr;
T * dev_c = nullptr;

cudaMalloc( (void**)&dev_a, size*sizeof(T) );
cudaMalloc( (void**)&dev_b, size*sizeof(T) );
cudaMalloc( (void**)&dev_c, size*sizeof(T) );

cudaMemcpy( dev_a, a, size*sizeof(T), cudaMemcpyHostToDevice );
cudaMemcpy( dev_b, b, size*sizeof(T), cudaMemcpyHostToDevice );

// 2 blocks ; (size+1)/2 threads per block 
addKernel<<< 2, (size + 1)/2 >>>(dev_a,dev_b,dev_c,size);

cudaDeviceSynchronize();

cudaMemcpy(c,dev_c,size*sizeof(T), cudaMemcpyDeviceToHost );
cudaFree(dev_c);
cudaFree(dev_b);
cudaFree(dev_a);
}

int main(void)
{
const int arraySize = 5;
const double a[arraySize] = {1,2,3,4,5};
const double b[arraySize] = {10,20,30,40,50};
double c[arraySize] = {0};

cudaAdd<double>(a,b,c,arraySize);

std::cout << "result = { " << c[0] << "; " << c[1] << "; "<< c[2] << "; "<< c[3] << "; "<< c[4] << " }" << std::endl;
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}
