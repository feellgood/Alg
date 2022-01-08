#ifndef GPU_COEFF_H
#define GPU_COEFF_H

namespace GPU
{
/**
\class v_coeff
container for pairs of indice nd double value to represent a coefficient of a sparse vector
*/
class v_coeff
{
public:
	/** no parameter constructor */
	inline __host__ __device__ v_coeff(): _i(0),_c(0.0) {}
	/** constructor */
	inline __host__ __device__ v_coeff(const size_t i,const double c):_i(i),_c(c) {}
	
	/** getter for the value of the coefficient */
	inline double getVal(void) const {return _c;} 
	
	/** setter for the value of the coefficient */
	inline void setVal(const double val) {_c = val;}

/** ref to coeff value */
inline double & valRef(void) {return _c;}	

	/** increment value with val */
	inline void inc(const double val) { _c += val;}
/**
lexicographic order
*/
	inline bool operator< (const v_coeff &c) const 
	{ return (this->_i < c._i); }

/**
two coeffs are equal when their indices are equals
*/
	inline bool operator== (const v_coeff &c) const
	{ return (this->_i == c._i); }

	/** index of the coefficient */
	size_t _i;

private:
/** value of the coefficient */	
	double _c;
}; //end class v_coeff

}

#endif
