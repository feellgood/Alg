#include "GPU_core.h"

int main(void)
{
thrust::host_vector<double> X(4);
X[0]=1.01;
X[1]=2.5;
X[2]=3.14;
X[3]=-1.23456789;

thrust::device_vector<double> dX = X; 
dX[3] += 1;

thrust::device_vector<double> dY(4);
dY[0]=1.0;
dY[1]=2.0;
dY[2]=-1.0;
dY[3]=0.0;

std::cout <<"X[3] on host : " << X[3] << std::endl;
std::cout <<"X[3] on device : " << dX[3] << std::endl;

std::cout << "X.Y = " << GPU::dot(dX,dY) << std::endl;
/*
GPU::scaled(dX,2.0,dY);
std::cout << "result = { " << dY[0] << "; " << dY[1] << "; "<< dY[2] << "; "<< dY[3] << "}" << std::endl;
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;

thrust::device_vector<double> dZ(4);

GPU::p_direct(dX,dY,dZ);
std::cout << "result = { " << dZ[0] << "; " << dZ[1] << "; "<< dZ[2] << "; "<< dZ[3] << "}" << std::endl;
*/
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;

return 0;
}
