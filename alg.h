#ifndef ALG_H
#define ALG_H

/** \file alg.h 
 * \brief set of class to handle sparse matrix operations for gradient conjugate algorithm
 * a sparse vector class
 * a read and a write sparse matrix class
 * various vector operations, scalar and direct products; ...
 * most of the member functions of the classes are using lambdas and C++11 algorithm and numeric
 * */

#include <iostream>
#include <cmath> // sqrt,fabs
#include <algorithm>
#include <numeric> // inner_product

#include "alg_coeff.h"
#include "alg_sparseVect.h"
#include "alg_sparseMat.h"
#include "alg_denseMat.h"
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

/** euclidian norm of vector X */
inline double norm(const std::vector<double> & X)
	{ return sqrt(fabs( p_scal(X,X) )); }

/** Y = A*X */
inline void mult(alg::r_sparseMat & A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size)
	{ for(size_t i=0;i<_size;i++) { Y[i]= alg::p_scal(A(i),X); } }
}

/** C = A*B*/
inline void mult(alg::denseMat const& A,alg::denseMat const& B,alg::denseMat & C)
{
const size_t nrowsA = A.nrows();
const size_t ncolsA = A.ncols();
const size_t ncolsB = B.ncols();
assert(ncolsA==B.nrows());

C.nrows()=nrowsA;
C.ncols()=ncolsB;
C.resize(nrowsA*ncolsB);
for (size_t i=0; i<nrowsA; i++) { 
    for (size_t j=0; j<ncolsB; j++) {
        C(i, j) = 0.0;
        for (size_t k=0; k<ncolsA; k++) { C(i,j) += A(i, k)*B(k, j); }
      } 
    }
}

/** conjugate gradient with diagonal preconditionner */
void cg_dir(r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, const std::vector<size_t>& ld, alg::iteration &iter);
} // end namespace alg

#endif //ALG_H

