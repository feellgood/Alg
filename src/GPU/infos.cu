#include <iostream>
#include <cuda_runtime.h>

void infos(void)
	{ cudaDeviceProp devProp;
	if (cudaSuccess == cudaGetDeviceProperties(&devProp,0))
	std::cout << "CUDA Device = " << devProp.major << "." << devProp.minor << " has "<< devProp.multiProcessorCount << " multi-processors" << std::endl;
	}
