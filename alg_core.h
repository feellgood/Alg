#ifndef ALG_CORE_H
#define ALG_CORE_H

/** \namespace alg
 * grab altogether many algebraic functions between various vectors and/or matrices
*/

#include <numeric> // inner_product
#include "alg_sparseMat.h"
#include "alg_denseMat.h"

namespace alg
{
/**  Y = alpha*X */
inline void scaled(const std::vector<double> & X, const double alpha, std::vector<double> & Y) 
	{ std::transform(X.begin(),X.end(),Y.begin(),[alpha](const double _x){ return alpha*_x; }); }

/**  Y *= alpha  */
inline void scaled( const double alpha, std::vector<double> & Y) 
	{ std::for_each(Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }

/** returns scalar product X.Y */
inline double dot(const std::vector<double> & X,const std::vector<double> & Y)
	{ return std::inner_product(X.begin(),X.end(),Y.begin(),0.0); }

/** direct product : Z = X⊗Y */
inline void p_direct(const std::vector<double> & X,const std::vector<double> & Y,std::vector<double> & Z)
	{ for(unsigned int i=0;i<Z.size();i++) Z[i]=X[i]*Y[i]; }

/** Y += X       */
inline void add(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::plus<double>()  ); }

/** Y -= X       */
inline void sub(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::minus<double>()  ); }

/** Y += alpha*X       */
inline void scaled_add(const std::vector<double> & X,const double alpha, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),[alpha] (const double _x,double _y) { return _x+(alpha*_y); }   ); }

/** euclidian norm of vector X */
inline double norm(const std::vector<double> & X)
	{ return sqrt(fabs( alg::dot(X,X) )); }

/** Y = A*X with sparseMat A */
inline void mult(alg::r_sparseMat & A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size)
	{ for(size_t i=0;i<_size;i++) Y[i]= A(i).dot(X); }
}

/** Y = A*X with denseMat A */
inline void mult(alg::denseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
const size_t ncolsA = A.ncols();
const size_t nrowsA = A.nrows();
assert(ncolsA==_size);

Y.resize(_size);
for (size_t i=0; i<nrowsA; i++) {
       double val(0);
       for (size_t j=0; j<_size; j++) val += A(i,j)*X[j];
       Y[i]=val;
       }
}

/** Y = trans(X)*A */
inline void transposed_mult(std::vector<double> const& X,alg::denseMat const& A,std::vector<double> &Y)
{
const size_t ncolsA = A.ncols();
Y.resize(ncolsA);
const size_t _size = X.size();

for (size_t j=0; j<ncolsA; j++) { 
       double val(0);
       for (size_t i=0; i<_size; i++) val += X[i]*A(i,j);
       Y[j]=val;
       }
}

/** Y = trans(A)*X with denseMat A */
inline void transposed_mult(alg::denseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
const size_t ncolsA = A.ncols();
const size_t nrowsA = A.nrows();
assert(nrowsA==_size);

Y.resize(ncolsA);
for (size_t i=0; i<ncolsA; i++) {
       double val(0);
       for (size_t j=0; j<nrowsA; j++) val += A(j, i)*X[j];
       Y[i]=val;
       }
}

/** C = trans(A)*B */
inline void transposed_mult(alg::denseMat const& A,alg::denseMat const& B,alg::denseMat & C)
{
const size_t nrowsA = A.nrows();
const size_t ncolsA = A.ncols();
const size_t ncolsB = B.ncols();
assert(nrowsA==B.nrows());

C.nrows()=ncolsA;
C.ncols()=ncolsB;
C.resize(ncolsA*ncolsB);
for (size_t i=0; i<ncolsA; i++) { 
    for (size_t j=0; j<ncolsB; j++) 
		{
	        double val(0);
	        for (size_t k=0; k<nrowsA; k++) { val += A(k, i)*B(k, j); }
		C(i, j) = val;      
		}
	}
}

/** C = A*B */
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
    for (size_t j=0; j<ncolsB; j++) 
	{
        double val(0);
        for (size_t k=0; k<ncolsA; k++) { val += A(i, k)*B(k, j); }
	C(i, j) = val;      
	}
     }
}

} // end namespace alg

#endif
