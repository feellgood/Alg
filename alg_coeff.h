#ifndef ALG_COEFF_H
#define ALG_COEFF_H

/** \file alg_coeff.h 
 * \brief set of class to handle sparse matrix and vector coefficients
 * two dedicated classes vector and matrix coefficients 
 * */

namespace alg
{
/**
\class v_coeff
container for pairs of indice nd double value to represent a coefficient of a sparse vector
*/
class v_coeff
{
public:
	/** constructor */
	inline v_coeff(const size_t i,const double c):_i(i),_c(c) {}
	
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


/**
\class m_coeff
container for a couple of indices and a double value to represent a coefficient of a sparse matrix
*/
class m_coeff
{
public:
	/** constructor */
	inline m_coeff(const size_t i,const size_t j,const double c):_i(i),_j(j),_c(c) {}
	
	/** getter for the value of the coefficient */
	inline double getVal(void) const {return _c;} 
/**
lexicographic order
*/
	inline bool operator< (const m_coeff &c) const 
	{ return ( (this->_i < c._i)||( (this->_i == c._i)&& (this->_j < c._j) ) ); }

/**
two coeffs are equal when their indices are equals
*/
	inline bool operator== (const m_coeff &c) const
	{ return ((this->_i == c._i)&&(this->_j == c._j)); }

	/** first index */
	size_t _i;

	/** second index */
	size_t _j;

private:
	/** value of the coefficient */	
	double _c;
}; //end class m_coeff 

} // end namespace alg


#endif //ALG_COEFF_H
