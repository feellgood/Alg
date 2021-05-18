#ifndef ALG_SPARSEMAT_H
#define ALG_SPARSEMAT_H

/** \file alg_sparseMat.h 
 \brief read and write sparse matrix
r_sparseMat : read sparse matrix
w_sparseMat : write sparse matrix : a std::vector of m_coeff = triplet of two indices and a value (i,j,value)
 */

#include <vector>

namespace alg
{

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

} // end namespace alg

#endif //ALG_SPARSEMAT_H
