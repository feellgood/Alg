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
void scaled( const double alpha, thrust::device_vector<double> & Y) 
	{
	auto op = [alpha] __device__ (double &_x){ _x *= alpha; };
	thrust::for_each(Y.begin(),Y.end(), op );
	}


/** return scalar product X.Y */
double dot(const thrust::device_vector<double> & X,const thrust::device_vector<double> & Y)
	{return thrust::inner_product(X.begin(),X.end(),Y.begin(),0.0); }

/** direct product : Z = X⊗Y */
void p_direct(const thrust::device_vector<double> & X,const thrust::device_vector<double> & Y,thrust::device_vector<double> & Z)
	{
	thrust::transform(X.begin(),X.end(),Y.begin(),Z.begin(),thrust::multiplies<double>()); 
	}

/** Y += X       */
void add(const thrust::device_vector<double> & X, thrust::device_vector<double> & Y)
	{ thrust::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),thrust::plus<double>() ); }

/** Y -= X       */
void sub(const thrust::device_vector<double> & X, thrust::device_vector<double> & Y)
	{ thrust::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),thrust::minus<double>() ); }

/** Y += alpha*X       */
void scaled_add(const thrust::device_vector<double> & X,const double alpha, thrust::device_vector<double> & Y)
	{
	auto op = [alpha] __device__ (const double _x,double _y) { return _x+(alpha*_y); };
	thrust::transform(Y.begin(),Y.end(),X.begin(),Y.begin(), op );
	}

/** euclidian norm of vector X */
double norm(const thrust::device_vector<double> & X)
	{ return sqrt(fabs( GPU::dot(X,X) )); }

} // end namespace GPU
	
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
*/
thrust::device_vector<double> dZ(4);

GPU::p_direct(dX,dY,dZ);
std::cout << "result = { " << dZ[0] << "; " << dZ[1] << "; "<< dZ[2] << "; "<< dZ[3] << "}" << std::endl;
std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;

return 0;
}
