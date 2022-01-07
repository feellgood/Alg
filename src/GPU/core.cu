#include<thrust/transform.h>
#include<thrust/for_each.h>
#include<thrust/inner_product.h>
#include<thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include<iostream>

namespace GPU
{

/**  Y = alpha*X */
void scaled(const thrust::device_vector<double> & X, const double alpha, thrust::device_vector<double> & Y) 
	{ 
	auto op = [alpha] __device__ (const double _x){ return alpha*_x; };
	
	thrust::transform(X.begin(),X.end(),Y.begin(), op ); 
	}

/**  Y *= alpha  */
/*
void scaled( const double alpha, thrust::device_vector<double> & Y) 
	{ thrust::for_each(Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }
*/

/** return scalar product X.Y */
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

GPU::scaled(dX,2.0,dY);
std::cout << "result = { " << dY[0] << "; " << dY[1] << "; "<< dY[2] << "; "<< dY[3] << "}" << std::endl;
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;

return 0;
}
