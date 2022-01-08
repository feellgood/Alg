#ifndef GPU_SPARSEVECT_H
#define GPU_SPARSEVECT_H

#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "GPU_coeff.h"

namespace GPU
{
/**
\class sparseVect
sparse vector : it is a container for v_coeff
*/
class sparseVect
{
public:
	/** dummy constructor */
	inline sparseVect() : sorted(true) {}

	/** constructor by initialization list */
	inline sparseVect(thrust::device_vector<GPU::v_coeff> _v) : sorted(false)
		{ x.assign(_v.begin(),_v.end()); } 

	/** inserter with value of a vector coefficient */
	inline void __host__ __device__ push_back(const size_t idx,const double c) { x.push_back(GPU::v_coeff(idx,c) ); sorted = false; }

	/** inserter with a vector coefficient */
	inline void __host__ __device__ push_back(GPU::v_coeff &coeff) { x.push_back(coeff); sorted = false; }

	/** sort the coefficients with lexicographic order */
	inline void sort() {thrust::sort(x.begin(),x.end()); sorted = true;}

	/** getter for sorted */
	inline bool isSorted(void) const {return sorted;}
	
	/** setter for sorted */
	inline void setSorted(bool s) {sorted = s;}
	
	/** getter for emptyness of the container of the coeffs */
	inline bool isEmpty(void) const {return x.empty();} 

	/** getter for the value of a coefficient of index idx, if several coeffs have the same index then it returns the value of the first occurence */
/*
	inline double getVal(size_t idx) const
		{
		double val(0);
		auto op = [idx] __device__ (GPU::v_coeff const& coeff){return (coeff._i == idx); };
		auto it = thrust::find_if(x.begin(),x.end(),op ); 
		if (it != x.end()) 
			{ 
			size_t i = thrust::distance(x.begin(),it); 
			GPU::v_coeff const& v_c = x[i];
			val = v_c.getVal(); 
			}		
		return val;		
		}
*/

	/** scalar product */
	inline double dot(const thrust::device_vector<double> & X) const
	{
	double val(0);
	const unsigned int X_dim = X.size();
//	auto op = [&X,X_dim,&val] __device__ (GPU::v_coeff coeff){ if(coeff._i < X_dim ) { val += coeff.getVal()*X.[coeff._i]; } };
//	thrust::for_each(x.begin(),x.end(), op);
		
	return val;
	}

private:

	/** coeffs container */
thrust::device_vector< GPU::v_coeff > x;

	/** if true the coeffs are sorted */
bool sorted; 

}; // end class sparseVect

/** scalar product of a sparse vector and a dense vector */
inline double dot(sparseVect const& X,const thrust::device_vector<double> & Y) { return X.dot(Y); }

}
#endif

