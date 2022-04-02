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

double xx[5] = {1.0f,2.0f,3.0f,4.0f,5.01f};
thrust::device_vector<double> dev_x(xx,xx+5);

/*
-1 * 1 + 1000.01 * 2 + 1 * 3 + 10 * 5.01 = 2052.12
*/
std::cout << "result1; result2 = " << toto.dot(dev_x) << "; " << GPU::dot(toto,dev_x) << std::endl;

double result = toto.getVal(1);
std::cout << "result = " << result  << std::endl;

//double y[5] = {1.0f,2.0f,0.0f,-1.0f,3.14f};


std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}

