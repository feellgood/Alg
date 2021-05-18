#ifndef ALG_H
#define ALG_H

/** \file alg.h 
 * \brief set of class to handle sparse matrix operations for gradient conjugate algorithm
 * a sparse vector class
 * a write sparse matrix class
 * a read sparse matrix class
 * various vector operations, scalar and direct products; ...
 * most of the member functions of the classes are using lambdas and C++11 algorithm and numeric
 * */

#include <iostream>
#include <cmath> // sqrt,fabs
#include <algorithm>
#include <numeric> // inner_product

#include "alg_coeff.h"
#include "alg_sparseVect.h"
#include "alg_iter.h"

/** \namespace alg
 * grab altogether sparse matrix, vector and dedicated functions
*/

namespace alg
{
/** 
Y = alpha*X 
*/
inline void scaled(const std::vector<double> & X, const double alpha, std::vector<double> & Y) 
	{ std::transform(X.begin(),X.end(),Y.begin(),[alpha](const double _x){ return alpha*_x; }); }

/** 
Y *= alpha 
*/
inline void scaled( const double alpha, std::vector<double> & Y) 
	{ std::for_each(Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }

/**
returns scalar product X.Y
*/
inline double p_scal(const std::vector<double> & X,const std::vector<double> & Y)
	{ return std::inner_product(X.begin(),X.end(),Y.begin(),0.0); }



/**
direct product : Z = X⊗Y
*/
inline void p_direct(const std::vector<double> & X,const std::vector<double> & Y,std::vector<double> & Z)
	{ for(unsigned int i=0;i<Z.size();i++) Z[i]=X[i]*Y[i]; }

/** Y += X       */
inline void inc(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::plus<double>()  ); }

/** Y -= X       */
inline void dec(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::minus<double>()  ); }


/** Y += alpha*X       */
inline void scaled_inc(const std::vector<double> & X,const double alpha, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),[alpha] (const double _x,double _y) { return _x+(alpha*_y); }   ); }


/**
euclidian norm of vector X
*/
inline double norm(const std::vector<double> & X)
	{ return sqrt(fabs( p_scal(X,X) )); }




/**
\class w_sparseMat
write sparse Matrix, it is a container for objects m_coeff. 
If some m_coeff have the same indices, they will be summed to build the real matrix coefficient using rebuild member function.(not implemented yet) 
*/
class w_sparseMat
{
	friend class r_sparseMat;

public:
	/** constructor */
	inline w_sparseMat(const size_t _N):N(_N) { sorted = false; collected = false; }
	
	/** inserter for a coefficient */
	inline void push_back(const m_coeff &co) { C.push_back(co); }
	
	/** inserter with direct values of a coefficient */
	inline void push_back(const size_t i,const size_t j, const double val) {C.push_back(alg::m_coeff(i,j,val));}
	
	/** getter for the number of lines */
	inline size_t getDim(void) const {return N;}

	/** sort the coefficients in lexicographic order and refresh collected and sorted booleans */
	inline void rebuild(void) 
		{ 
		std::sort(C.begin(),C.end());
	       	sorted = true;	
		std::vector<m_coeff>::iterator it = std::adjacent_find(C.begin(),C.end() );		
		if (it == C.end()) { collected = true; } 
		else { collected = false; }		
		}

	/** getter for sorted */
	inline bool isSorted(void) const {return sorted;}
	
	/** getter for collected */
	inline bool isCollected(void) const {return collected;}

	/** printing function */
	inline void print(std::ostream & flux) const
	{ flux<<'{'; std::for_each(C.begin(),C.end(), [&flux](const m_coeff &c){ flux << '{' << c._i << ','<< c._j << ':' << c.getVal() <<'}';}); flux<<"}\n"; }


private:
/** dimension of sparse matrix, N is the number of lines */
	size_t N;

/** if sorted == true, coeffs have been sorted in lexicographic order */
	bool sorted;

	/** if collected == true, coefficients have been regrouped in lexicographic order,and redundant coeffs summed (if any) */
	bool collected;

	/**
container for the sparse matrix coefficient, C.size() might be different from N, since a coefficient with indices (i,j) might be push_backed several times
*/
	std::vector<alg::m_coeff> C;
}; // end class w_sparseMat


/** operator<< for w_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, w_sparseMat const& m) {m.print(flux); return flux;}

/** operator<< for sparseVect */
inline std::ostream & operator<<(std::ostream & flux, sparseVect const& v) {v.print(flux); return flux;}

/** scalar product of a sparse vector and a dense vector */
inline double p_scal(sparseVect const& X,const std::vector<double> & Y)
	{ return X.p_scal(Y); }

	/** \class r_sparseMat
read sparse matrix	 
	The constructor is buiding from a write sparse matrix the data to access efficiently the coefficients values
       	*/
class r_sparseMat
{
public:
	/** constructor */
	inline r_sparseMat(w_sparseMat &A):N(A.getDim())
		{
		m.resize(N);// N is the number of lines
		if (!A.isSorted()) { A.rebuild(); }
		
		if (!A.C.empty())
			{
			for(std::vector<m_coeff>::iterator it = A.C.begin(); it != A.C.end() ; ++it)
				{ if (it->_i < N) { m[it->_i].push_back(it->_j,it->getVal());} }	
			
			std::for_each(m.begin(),m.end(),[](sparseVect & _v) {_v.collect();} );
			}
		}
	/** printing function */
	inline void print(void) { std::for_each(m.begin(),m.end(),[](sparseVect const& _v) {std::cout << _v;} ); }

	/** getter for the number of lines */
	inline size_t getDim(void) const {return N;}

	/** getter for an innner sparse vector */
	inline alg::sparseVect & operator() (const size_t & i) {return m[i];}

	/** getter for a coefficient value */
	inline double operator() (const size_t &i, const size_t &j) { return m[i].getVal(j); }

private:
/** dimension of the sparse matrix (nb of lines) */
	const size_t N;

	/** coefficient container */
	std::vector<sparseVect> m;
}; // end class r_sparseMat

/** Y = A*X */
inline void mult(alg::r_sparseMat & A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size)
	{ for(size_t i=0;i<_size;i++) { Y[i]= alg::p_scal(A(i),X); } }
}

/** conjugate gradient with diagonal preconditionner */
void cg_dir(r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, const std::vector<size_t>& ld, alg::iteration &iter);
} // end namespace alg

#endif //ALG_H

