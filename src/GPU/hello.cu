#include<stdio.h>
#include <iostream>

__global__ void hello_kernel(void)
{ 
printf("Hello World! (from GPU)\n");// we have to speak C here (GPU)
}

int main(void)
{
std::cout << "Hello World! (from host)\n"; // here (host) we can use C++
hello_kernel<<<1,1>>>();

//wait for device to finish to see the result (here the 'hello world' from the GPU)
cudaDeviceSynchronize();
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}
