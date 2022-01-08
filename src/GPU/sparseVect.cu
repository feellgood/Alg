#include "GPU_core.h"

#include "GPU_sparseVect.h"

#include<iostream>

int main(void)
{
GPU::sparseVect toto;
toto.push_back(0,1.23);

std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}

