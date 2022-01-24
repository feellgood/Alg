#include "GPU_core.h"

#include "GPU_sparseVect.h"

#include<iostream>
#include<iomanip>

int main(void)
{
GPU::sparseVect toto;

toto.push_back(0, -1.0f);

toto.push_back(2, 1.0f);
toto.push_back(4, 10.0f);

toto.push_back(1, 1000.01f);

//toto.sort();

double x[5] = {1.0f,2.0f,3.0f,4.0f,5.01f};
thrust::device_vector<double> dev_x(x,x+5);

/*
-1 * 1 + 1000.01 * 2 + 1 * 3 + 10 * 5.01 = 2052.12
*/
std::cout << "result = " << toto.dot(dev_x) << std::endl;

std::cout << "result = " << GPU::dot(toto,dev_x) << std::endl;
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}

