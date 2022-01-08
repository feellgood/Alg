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


