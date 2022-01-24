#ifndef GPU_SPARSEVECT_H
#define GPU_SPARSEVECT_H

#include <thrust/device_vector.h>
#include <thrust/sort.h>

//#include "GPU_coeff.h"

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
//	inline sparseVect(thrust::device_vector<GPU::v_coeff> _v) : sorted(false) { x.assign(_v.begin(),_v.end()); } 

	/** inserter with value of a vector coefficient */
	inline void __host__ __device__ push_back(const size_t _idx,const double c) { idx.push_back(_idx); x.push_back(c); sorted = false; }

	/** sort the coefficients with order greater_than on indices ; both idx and x are sorted */
	inline void sort() {thrust::sort_by_key(idx.begin(),idx.end(),x.begin() ); sorted = true;}

	/** getter for sorted */
	inline bool isSorted(void) const {return sorted;}
	
	/** setter for sorted */
	inline void setSorted(bool s) {sorted = s;}
	
	/** getter for emptyness of the container of the coeffs */
	inline bool isEmpty(void) const {return x.empty();} 

	/** getter for the value of a coefficient of index idx, if several coeffs have the same index then it returns the value of the first occurence */

	double getVal(const size_t _idx) const
		{
		double val(0);
		auto result = thrust::find(idx.begin(),idx.end(), _idx ); 
		
		if (result != idx.end()) 
			{
			size_t _i = thrust::distance(idx.begin(),result);  
			val = x[ _i ];
			}		

		return val;		
		}


	/** scalar product */
	double dot(const thrust::device_vector<double> & X) const
	{
	return thrust::inner_product( thrust::make_permutation_iterator(X.begin(),idx.begin()),
	thrust::make_permutation_iterator(X.begin(),idx.end()), x.begin(), 0.0f ); 
	// carefull here : init  must be 0.0f if 0 there is a cast to int that propagates, and result is wrong
	}

private:

	/** coeffs container */
thrust::device_vector< size_t> idx;
thrust::device_vector< double > x;

	/** if true the coeffs are sorted */
bool sorted; 

}; // end class sparseVect

/** scalar product of a sparse vector and a dense vector */
double dot(sparseVect const& X,const thrust::device_vector<double> & Y) { return X.dot(Y); }

}
#endif

