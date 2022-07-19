#ifndef ALG_SPARSEVECT_H
#define ALG_SPARSEVECT_H

/** \file alg_sparseVect.h 
 * \brief sparse vector
a sparse vector is a derived class from SparseVect abstract class. It is a std::vector of v_coeff, which is a couple composed of an index and a value
to populate with coefficients the sparse vector, use push_back method
If several v_coeff have the same index, then they are summed when calling the method collect()
It is possible to kill v_coeffs with index idx, calling kill(idx), it will erase all v_ceffs with index idx
It is also possible to erase v_coeffs with value zero calling kill_zero, it may be usefull to spare some memory.  
 * */

#include <execution>
#include <vector>
#include <iostream> 

#include "alg_coeff.h"

namespace alg
{
/** 
\class SparseVect
Abstract sparse vector class, push_back to insert in a container the coefficients, dot to compute scalar product.
*/
class SparseVect
{
public:
/** push_back method to insert a coefficient in the sparse vector container. If the coefficient does not exist it has to be created, if it does the coeff should be either += added or included in a container to treat it later, with sort() and collect() methods if needed  */
	virtual void push_back(const size_t idx,const double c) = 0;

	/** return the value of the coefficient of index idx, if idx is not unique returns the first occurence */
	virtual double getVal(const size_t idx) const = 0;

/** scalar product */
	virtual double dot(std::vector<double> const& X) const = 0;
};



/**
\class sparseVect
sparse vector : container for v_coeff
*/

//class sparseVect
class sparseVect : public SparseVect
{
public:
	/** dummy constructor */
	inline sparseVect() : sorted(true),collected(false) {}

	/** constructor by initialization list */
	inline sparseVect(std::vector<alg::v_coeff> _v) : sorted(false),collected(false) { x.assign(_v.begin(),_v.end()); } 

	/** inserter with value of a vector coefficient */
	inline void push_back(const size_t idx,const double c) { x.push_back(alg::v_coeff(idx,c) ); sorted = false; collected=false; }

	/** inserter with a vector coefficient */
	inline void push_back(alg::v_coeff &coeff) { x.push_back(coeff); sorted = false; collected = false; }

	/** sort the coefficients with lexicographic order */
	inline void sort() {std::sort(std::execution::par_unseq,x.begin(),x.end()); sorted = true;}

	/** erase all coefficients with index idx */
	inline void kill(const size_t idx) 
		{x.erase(std::remove_if(x.begin(),x.end(),[this,&idx](alg::v_coeff &c) { return (c._i == idx); } ),x.end());}

	/** erase all coefficients with value zero */
	inline void kill_zero(void) 
		{x.erase(std::remove_if(x.begin(),x.end(),[this](alg::v_coeff &c) { return (c.getVal() == 0); }),x.end() );}

/** return true if the coefficient exists */
	inline bool exist(const size_t &idx) const
		{ 
		auto it = std::find_if(x.begin(),x.end(),[this,&idx](alg::v_coeff coeff){return (coeff._i == idx); } ); 
		return (it != x.end());
		}

	/** collect method is sorting all v_coeffs, eventually with redundant indices, and is summing coeffs with same indices. It removes the coeffs that have been summed. */
	inline void collect(void)
		{
		if (!sorted) sort();
		std::vector<v_coeff>::iterator it = std::adjacent_find(std::execution::par,x.begin(),x.end());
		while (it != x.end())
			{
			alg::v_coeff c(it->_i,0);			
			while((it != x.end()) && (c == *it)  )			
				{ c.inc( it->getVal() ); ++it; }			
			kill(c._i);			
			push_back(c);			
			it = std::adjacent_find(x.begin(),x.end());			
			}		
		collected = true;		
		}

	/** getter for sorted */
	inline bool isSorted(void) const {return sorted;}
	
	/** getter for emptyness of the container of the coeffs */
	inline bool isEmpty(void) const {return x.empty();} 

	/** getter for collected */
	inline bool isCollected(void) const {return collected;}

	/** getter for the value of a coefficient of index idx, if several coeffs have the same index then it returns the value of the first occurence */
	double getVal(const size_t idx) const
		{
		double val(0);
		auto pred = [&idx](alg::v_coeff const& coeff){return (coeff._i == idx); };
		auto it = std::find_if(x.begin(),x.end(), pred ); 
		if (it != x.end()) val = it->getVal();		
		return val;		
		}

	/** return a reference to the value of the coefficient idx, we should never do that ... */
	/*
	inline double & getValRef(size_t idx)
		{
		auto it = std::find_if(x.begin(),x.end(),[&idx](alg::v_coeff & coeff){return (coeff._i == idx); } ); 
		return it->valRef();// carefull might be out of bounds when it == x.end()	
		}
	*/

	/** setter for the value of a coefficient of index idx, all coeffs must have a unique idx, call collect() method before if needed */
	/*
	inline void setVal(const size_t idx,const double val)
		{
		if (collected)
			{
			auto it = std::find_if(x.begin(),x.end(),[&idx](alg::v_coeff & coeff){return (coeff._i == idx); } );
			if (it != x.end()) it->setVal(val);
			}
		}
	*/

	/** scalar product, callable whatever collected and sorted booleans are */
	double dot(std::vector<double> const& X) const
	{
	double val(0);
	
	auto op = [&X,&val](alg::v_coeff const& coeff){ val += coeff.getVal()*X.at(coeff._i);  }; // could be replaced by accumulate
	std::for_each(x.begin(),x.end(), op );
		
	return val;
	}

	/** printing function */
	inline void print(std::ostream & flux) const
	{ flux<<'{'; std::for_each(x.begin(),x.end(), [&flux](const v_coeff &c){ flux << '{' << c._i << ':' << c.getVal() <<'}';}); flux<<"}\n"; }

	/** iterator : begin */
	inline std::vector<alg::v_coeff>::iterator begin() { return x.begin(); }

	/**iterator : end */
	inline std::vector<alg::v_coeff>::iterator end() { return x.end(); }

	/** const_iterator : begin */
	inline std::vector<alg::v_coeff>::const_iterator begin() const { return x.begin(); } 

	/** const_iterator : end */
	inline std::vector<alg::v_coeff>::const_iterator end() const { return x.end(); } 

	/** const_iterator : cbegin */
	inline std::vector<alg::v_coeff>::const_iterator cbegin() const { return x.cbegin(); } 

	/** const_iterator : cend */
	inline std::vector<alg::v_coeff>::const_iterator cend() const { return x.cend(); } 

private:

	/** coeffs container */
std::vector< alg::v_coeff > x;

	/** if true the coeffs are sorted */
bool sorted; 

/** if true the coeffs have been collected */
bool collected;
}; // end class sparseVect

/** scalar product of a sparse vector and a dense vector */
inline double dot(sparseVect const& X,const std::vector<double> & Y) { return X.dot(Y); }

/** operator<< for sparseVect */
inline std::ostream & operator<<(std::ostream & flux, sparseVect const& v) {v.print(flux); return flux;}

} // end namespace alg

#endif //ALG_SPARSEVECT_H
