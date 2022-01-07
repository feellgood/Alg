#include<thrust/inner_product.h>
#include<thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include<iostream>

namespace GPU
{
double dot(const thrust::device_vector<double> & X,const thrust::device_vector<double> & Y)
	{return thrust::inner_product(X.begin(),X.end(),Y.begin(),0.0); }
}
	
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

return 0;
}
