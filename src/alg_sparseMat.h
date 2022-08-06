#ifndef ALG_SPARSEMAT_H
#define ALG_SPARSEMAT_H

/** \file alg_sparseMat.h 
 \brief read and write sparse matrix
r_sparseMat : read sparse matrix : it is buit calling the constructor with a w_sparseMat as argument
w_sparseMat : write sparse matrix : a std::vector of m_coeff = triplet of two indices and a value (i,j,value)
 */

#include <vector>

#include "alg_sparseVect.h"
#include "alg_utils.h"


namespace alg
{

/**
\class SparseMat
class T must have print() method and dot method
*/

template <typename T,typename T_elem>
class SparseMat
{
public:
	/** insert the coefficient of index (i,j) */
	virtual void push_back(const size_t i, const size_t j, const T_elem c) = 0;

	/** getter for a coefficient value */
	virtual T_elem operator() (const size_t i, const size_t j) const = 0;//{ return M[i].getVal(j); }

	/** returns dimension N*/
	size_t getDim(void) const {return N;}


	/** printing function */
	void print(void) const
		{std::for_each(M.begin(),M.end(),[](T const& _v) {_v.print();} );}
		//{ std::for_each(M.begin(),M.end(),[](T const& _v) {std::cout << _v;} ); }

	/** printing function */
	void print(std::ostream & flux) const
		{ std::for_each(M.begin(),M.end(),[&flux](T const& _v) {_v.print(flux);} ); }
		
	/** matrix vector multiplication : Y = M * X */
	inline void mult(std::vector<T_elem> const& X,std::vector<T_elem> &Y) const
		{ std::transform(std::execution::par_unseq,M.begin(),M.end(),Y.begin(),[&X] (T const& _v) { return _v.dot(X); } ); }

protected:
	/** number of lines */
	size_t N;
	
	/** matrix coefficients */
	std::vector<T> M;
};



/**
\class w_sparseMat
write sparse Matrix, it is a container for objects m_coeff. 
If some m_coeff have the same indices, they will be summed to build the real matrix coefficient using rebuild member function.
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
	
	/** sort the coefficients in lexicographic order and refresh collected and sorted booleans */
	inline void rebuild(void) 
		{ 
		std::sort(std::execution::par_unseq,C.begin(),C.end());
	       	sorted = true;	
		auto it = std::adjacent_find(C.begin(),C.end() ); //std::vector<m_coeff>::iterator it = std::adjacent_find(M.begin(),M.end() );		
		if (it == C.end()) { collected = true; } 
		else { collected = false; }		
		}

/** returns dimension N*/
	size_t getDim(void) const {return N;}

	/** setter for sorted */
	inline void setSorted(bool b) {sorted = b;}

	/** setter for collected */
	inline void setCollected(bool b) {collected = b;}

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


/** \class sparseMat
read sparse matrix	 
	The constructor is buiding from a write sparse matrix the data to access efficiently the coefficients values
*/
class r_sparseMat : public SparseMat<sparseVect,double>
{
public:
	/** constructor */
	r_sparseMat(const size_t dim) { N = dim; M.resize(dim); }// dim is the number of lines

	/** constructor : build coefficients from write sparse matrix */
	r_sparseMat(w_sparseMat &A)
		{
		N = A.getDim();
		M.resize(N);// N is the number of lines
		
		if (!A.C.empty())
			{
			if (!A.isSorted()) { A.rebuild(); }
			for(auto it = A.C.begin(); it != A.C.end() ; ++it)
				{ if (it->_i < N) M[it->_i].push_back(it->_j,it->getVal()); }
			
			collect();
			}
		}

/** inserter with direct values of a coefficient */
	inline void push_back(const size_t i,const size_t j, const double val) {M[i].push_back(j,val);}

/** return true if the coefficient exists */
	inline bool exist(const size_t &i, const size_t &j) const { return ( (i<N)&&(M[i].exist(j)) ); }

	/** getter for an innner sparse vector */
	inline alg::sparseVect & operator() (const size_t & i) {return M[i];}

	/** getter for a coefficient value */
	inline double operator() (const size_t i, const size_t j) const { return M[i].getVal(j); }

	/** call collect method for all sparse vectors  */
	inline void collect(void) { std::for_each(std::execution::par_unseq,M.begin(),M.end(),[](sparseVect & _v) {_v.collect();} ); }


/** matrix vector multiplication : Y = this.X according mask b */
	inline void maskedMult(std::vector<bool> const &b, std::vector<double> const& X, std::vector<double> &Y) const
		{
		std::transform(std::execution::par_unseq,M.begin(),M.end(),b.begin(),Y.begin(),
		[&X] (alg::sparseVect const&_v,const bool _b) { return( _b ? _v.dot(X) : 0.0 ); } 
		); 
		}


/** build diagonal preconditionner */
void buildDiagPrecond(std::vector<double> &DP) const
	{ for(size_t i=0;i<N;i++) { DP[i] = 1.0/(M[i].getVal(i)); } }

/** build diagonal preconditionner respecting mask b */
void buildDiagPrecond(std::vector<bool> const& b,std::vector<double> &DP) const
	{ for(size_t i=0;i<N;i++) { DP[i] = (b[i] ? 1.0/(M[i].getVal(i)) : 0.0 ); } }

}; // end class r_sparseMat

/** operator<< for r_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, r_sparseMat const& m) {m.print(flux); return flux;}

/** operator<< for w_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, w_sparseMat const& m) {m.print(flux); return flux;}

} // end namespace alg

#endif //ALG_SPARSEMAT_H
