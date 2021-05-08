#ifndef ALG_H
#define ALG_H

#include <iostream>
#include <cmath> // sqrt,fabs
#include <algorithm>
#include <numeric> // inner_product

#include "alg_iter.h"

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
norme de X
*/
inline double norm(const std::vector<double> & X)
	{ return sqrt(fabs( p_scal(X,X) )); }


/**
\class v_coeff
container for pairs of indice nd double value to represent a coefficient of a sparse vector
*/
class v_coeff
{
public:
	inline v_coeff(const size_t i,const double c):_i(i),_c(c) {}
	inline double getVal(void) const {return _c;} 
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

	size_t _i;

private:	
	double _c;
}; //end class v_coeff


/**
\class m_coeff
container for a couple of indices and a double value to represent a coefficient of a sparse matrix
*/
class m_coeff
{
public:
	inline m_coeff(const size_t i,const size_t j,const double c):_i(i),_j(j),_c(c) {}
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

	size_t _i;
	size_t _j;

private:	
	double _c;
}; //end class m_coeff 

/**
\class w_sparseMat
write sparse Matrix, it is a container for objects m_coeff. 
If some m_coeff have the same indices, they will be summed to build the real matrix coefficient using rebuild member function.(not implemented yet) 
*/
class w_sparseMat
{
	friend class r_sparseMat;

public:
	inline w_sparseMat(const size_t _N):N(_N) { sorted = false; collected = false; }
	inline void push_back(const m_coeff &co) { C.push_back(co); }
	inline void push_back(const size_t i,const size_t j, const double val) {C.push_back(alg::m_coeff(i,j,val));}
	inline size_t getDim(void) const {return N;}

	inline void rebuild(void) 
		{ 
		std::sort(C.begin(),C.end());
	       	sorted = true;	
		std::vector<m_coeff>::iterator it = std::adjacent_find(C.begin(),C.end() );		
		if (it == C.end()) { collected = true; } 
		else { collected = false; }		
		}

	inline bool isSorted(void) const {return sorted;}
	inline bool isCollected(void) const {return collected;}

	inline void print(std::ostream & flux) const
	{ flux<<'{'; std::for_each(C.begin(),C.end(), [&flux](const m_coeff &c){ flux << '{' << c._i << ','<< c._j << ':' << c.getVal() <<'}';}); flux<<"}\n"; }


private:
/** dimension of sparse matrix, N is the number of lines */
	size_t N;
	bool sorted;
	bool collected;
/**
container for the sparse matrix coefficient, C.size() might be different from N, since a coefficient with indices (i,j) might be push_backed several times
*/
	std::vector<alg::m_coeff> C;
}; // end class w_sparseMat

/**
\class sparseVect
sparse vector : it is a container for v_coeff
*/
class sparseVect
{
public:
	inline sparseVect() {sorted = true;collected = false;}

	inline void push_back(const size_t idx,const double c) { x.push_back(alg::v_coeff(idx,c) ); sorted = false; collected=false; }

	inline void push_back(alg::v_coeff &coeff) { x.push_back(coeff); sorted = false; collected = false; }

	inline void sort() {std::sort(x.begin(),x.end()); sorted = true;}

	inline void kill(const size_t idx) 
		{x.erase(std::remove_if(x.begin(),x.end(),[this,&idx](alg::v_coeff &c) { return (c._i == idx); } ),x.end());}

	inline void kill_zero(void) 
		{x.erase(std::remove_if(x.begin(),x.end(),[this](alg::v_coeff &c) { return (c.getVal() == 0); }),x.end() );}

	/** collect method is sorting all v_coeffs, eventually with redundant indices, and is summing coeffs with same indices. It removes the coeffs that have been summed. */
	inline void collect(void)
		{
		if (!sorted) sort();
		std::vector<v_coeff>::iterator it = std::adjacent_find(x.begin(),x.end());
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

	inline bool isSorted(void) const {return sorted;}
	inline bool isEmpty(void) const {return x.empty();} 
	inline bool isCollected(void) const {return collected;}

	inline double getVal(size_t idx) const
		{
		double val(0);
		auto it = std::find_if(x.begin(),x.end(),[this,&idx](alg::v_coeff coeff){return (coeff._i == idx); } ); 
		if (it != x.end()) val = it->getVal();		
		return val;		
		}

	inline double p_scal(const std::vector<double> & X) const
	{
	double val(0);
	if (!isCollected()) {std::cout << "warning : cannot p_scal on an uncollected sparseVect" << std::endl;exit(1);}
	else
		{const unsigned int X_dim = X.size();
		for(auto it=x.begin();it!=x.end();++it)
			{ if(it->_i < X_dim ) { val += it->getVal()*X[it->_i]; } }
		}	
	return val;
	}

	inline void print(std::ostream & flux) const
	{ flux<<'{'; std::for_each(x.begin(),x.end(), [&flux](const v_coeff &c){ flux << '{' << c._i << ':' << c.getVal() <<'}';}); flux<<"}\n"; }

private:
std::vector< alg::v_coeff > x;
bool sorted; 
bool collected;
}; // end class sparseVect

/** operator<< for w_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, w_sparseMat const& m) {m.print(flux); return flux;}

/** operator<< for sparseVect */
inline std::ostream & operator<<(std::ostream & flux, sparseVect const& v) {v.print(flux); return flux;}

/** scalar product of a sparse vector and a dense vector */
inline double p_scal(sparseVect const& X,const std::vector<double> & Y)
	{ return X.p_scal(Y); }

class r_sparseMat
{
public:
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
	inline void print(void) { std::for_each(m.begin(),m.end(),[](sparseVect const& _v) {std::cout << _v;} ); }

	inline size_t getDim(void) const {return N;}

	inline alg::sparseVect & operator() (const size_t & i) {return m[i];}

	inline double operator() (const size_t &i, const size_t &j) { return m[i].getVal(j); }

private:
/** dimension of the square sparse matrix */
	const size_t N;
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
}

#endif //ALG_H

